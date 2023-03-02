# _ScreenModule.py
__module_name__ = "_ScreenModule.py"
__author__ = ", ".join(["Michael E. Vinyard", "Jayoung Kim Ryu"])
__email__ = ", ".join(["vinyard@g.harvard.edu", "jayoung_ryu@g.harvard.edu"])

import copy
import warnings
from typing import Union

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData

from .._arithmetic._funcs._log_fold_change import _log_fold_change
from .._normalization._funcs._read_count_norm import _log_normalize_read_count
from .._readwrite._funcs._read_screen_from_PoolQ import _read_screen_from_PoolQ
from .._readwrite._funcs._write_screen_to_csv import _write_screen_to_csv
from .._readwrite._funcs._write_screen_to_excel import _write_screen_to_excel
from .._utilities._funcs._update_dict import _update_dict
from ._supporting_functions._guides._GuideAnnotationModule import _annotate_sgRNAs
from ._supporting_functions._print_screen_object import _print_screen_object


class _Screen(AnnData):
    def __init__(self, X=None, guides=None, condit=None, *args, **kwargs):
        super().__init__(X, dtype=X.dtype, obs=guides, var=condit, *args, **kwargs)

    @property
    def guides(self):
        return self.obs

    @property
    def condit(self):
        return self.condit

    def __repr__(self) -> str:
        return _print_screen_object(self)[2]

    def __add__(self, other):
        if all(self.guides.index == other.guides.index) and all(
            self.condit.index == other.condit.index
        ):
            return _Screen(self.X + other.X, self.guides.copy(), self.condit.copy())
        else:
            raise ValueError("Guides/sample description mismatch")

    def read_PoolQ(self, path, metadata=False, merge_metadata_on="Condition"):
        """Read poolQ."""
        self._PoolQ_outpath = path
        self._PoolQScreenDict = _read_screen_from_PoolQ(self._PoolQ_outpath)

        for key, value in _update_dict(self._PoolQScreenDict).items():
            self.__setattr__(key, value)

        if metadata:
            self.condit = self.condit.merge(pd.read_csv(metadata), on=merge_metadata_on)

        _print_screen_object(self)

    def annotate_guides(
        self, genes, chrom, start, stop, annotations, DirectPairDict, ref_seq_path
    ):

        """
        Annotate sgRNA table.

        """
        self.guides = _annotate_sgRNAs(
            self.guides,
            genes,
            chrom,
            start,
            stop,
            annotations,
            DirectPairDict,
            ref_seq_path,
        )

    def log_norm(self, output_layer="lognorm_counts", read_count_layer=None):
        if read_count_layer is None:
            self.layers[output_layer] = _log_normalize_read_count(self.X)
        else:
            output_layer = f"lognorm_{read_count_layer}"
            self.layers[output_layer] = _log_normalize_read_count(
                self.layers[read_count_layer]
            )

    # TBD: mask ones with too low raw counts.
    def log_fold_change(
        self,
        cond1,
        cond2,
        lognorm_counts_key="lognorm_counts",
        name=False,
        out_guides_suffix="lfc",
        return_result=False,
    ):

        """
        General module to calculate LFC across experimental conditions.
        """
        if "lognorm" not in lognorm_counts_key:
            warnings.warn(
                "The layer specified must be log-normalized values using screen.log_norm()."
            )

        if lognorm_counts_key not in self.layers.keys():
            raise ValueError(
                "Specified normalized count isn't in your layer. First run screen.log_norm()."
            )

        cond1_idx = np.where(cond1 == self.condit.index)[0]
        cond2_idx = np.where(cond2 == self.condit.index)[0]
        if len(cond1_idx) != 1 or len(cond2_idx) != 1:
            if len(cond1_idx) == 0:
                print(f"No condition named {cond1} in Screen object.")
            else:
                print(f"Duplicate condition name {cond1} in Screen object")
            if len(cond2_idx) == 0:
                print(f"No condition named {cond2} in Screen object.")
            else:
                print(f"Duplicate condition name {cond2} in Screen object")
            raise ValueError("")

        try:
            lfc = _log_fold_change(
                self.layers[lognorm_counts_key], cond1_idx, cond2_idx
            )
            if return_result:
                return lfc
            else:
                self.guides[f"{cond1}_{cond2}.{out_guides_suffix}"] = lfc
        except:  # TBD: what error?
            print("Calculating LFC against two previously calculated LFC values...")
            dlfc = _log_fold_change(self.guides, cond1, cond2)

            if not name:
                name = (
                    f'{cond1.strip(".lfc")}_{cond2.strip(".lfc")}.d{out_guides_suffix}'
                )

            if return_result:
                return dlfc
            else:
                self.guides[name] = dlfc

    def log_fold_change_reps(
        self,
        cond1,
        cond2,
        lognorm_counts_key="lognorm_counts",
        rep_condit="replicate",
        compare_condit="sort",
        out_guides_suffix="lfc",
        keep_result=False,
    ):

        if rep_condit not in self.condit.columns:
            raise ValueError(f"{rep_condit} not in condit features")
        if compare_condit not in self.condit.columns:
            raise ValueError(f"{compare_condit} not in condit features")

        lfcs = []
        for rep in self.condit[rep_condit].unique():
            cond1_idx = np.where(
                (self.condit[rep_condit] == rep)
                & (self.condit[compare_condit] == cond1)
            )[0]
            cond2_idx = np.where(
                (self.condit[rep_condit] == rep)
                & (self.condit[compare_condit] == cond2)
            )[0]

            if len(cond1_idx) != 1 or len(cond2_idx) != 1:
                raise ValueError(
                    "Conditions are not unique for each replicates to be aggregated."
                )

            lfcs.append(
                self.log_fold_change(
                    self.condit.index[cond1_idx].tolist()[0],
                    self.condit.index[cond2_idx].tolist()[0],
                    lognorm_counts_key=lognorm_counts_key,
                    return_result=True,
                )
            )

        lfcs_array = np.concatenate(lfcs, axis=1)
        lfcs_df_columns = [
            f"{s}.{cond1}_{cond2}.{out_guides_suffix}"
            for s in self.condit[rep_condit].unique()
        ]
        lfcs_df = pd.DataFrame(
            lfcs_array, index=self.guides.index, columns=lfcs_df_columns
        )

        if keep_result:
            self.guides[lfcs_df_columns] = lfcs_df
        return lfcs_df

    # TODO: add guides metadata on how aggregates are calcualted?
    def log_fold_change_aggregate(
        self,
        cond1,
        cond2,
        lognorm_counts_key="lognorm_counts",
        aggregate_condit="replicate",
        compare_condit="sort",
        out_guides_suffix="lfc",
        aggregate_fn="median",
        name=None,
        return_result=False,
        keep_per_replicate=False,
    ):

        lfcs_df = self.log_fold_change_reps(
            cond1,
            cond2,
            lognorm_counts_key=lognorm_counts_key,
            rep_condit=aggregate_condit,
            compare_condit=compare_condit,
            out_guides_suffix=out_guides_suffix,
            keep_result=keep_per_replicate,
        )

        if aggregate_fn == "mean":
            lfcs_agg = lfcs_df.apply(np.mean, axis=1)
        elif aggregate_fn == "median":
            lfcs_agg = lfcs_df.apply(np.median, axis=1)
        elif aggregate_fn == "sd":
            lfcs_agg = lfcs_df.apply(np.std, axis=1)
        else:
            raise ValueError(
                "Only 'mean', 'median', and 'sd' are supported for aggregating LFCs."
            )

        if return_result:
            return lfcs_agg
        if name is None:
            self.guides[
                f"{cond1}_{cond2}.{out_guides_suffix}.{aggregate_fn}"
            ] = lfcs_agg
        else:
            self.guides[name] = lfcs_agg

    def fold_change(
        self,
        cond1: str,
        cond2: str,
        lognorm_counts_key: str = "lognorm_counts",
        return_result: bool = False,
    ) -> Union[None, pd.Series]:
        """Calculate log fold change (cond1/cond2) of normalized guide abundances."""

        log_fold_change = _log_fold_change(
            self.layers[lognorm_counts_key], cond1, cond2
        )
        if return_result:
            return log_fold_change
        self.guides[f"{cond1}_{cond2}.fc"] = log_fold_change

    def to_Excel(
        self,
        workbook_path: str = "CRISPR_screen.workbook.xlsx",
        index: bool = False,
        silent: bool = False,
        include_uns: bool = False,
    ) -> None:
        """Write components of Screen class to an Excel workbook.

        Args
        ---
        workbook_path: Prevent printing outpaths / details of the created workbook.
        index: If True, include an index in the workbook sheet.
        silent: If True, prevent printing outpaths / details of the created workbook.

        Notes:
        ------
        (1) Will likely need to be updated once we fully transition over to AnnData-like class.
        """

        _write_screen_to_excel(
            self,
            workbook_path,
            index,
            silent,
            include_uns,
        )

    def to_csv(self, out_path="CRISPR_screen"):

        """

        Write .csv files for each part of the screen. will eventually be replaced by something more native to AnnData.

        """

        _write_screen_to_csv(self, out_path)

    def to_mageck_input(
        self,
        out_path=None,
        count_layer=None,
        sgrna_column=None,
        target_column="target_id",
        sample_prefix="",
    ):
        """Formats screen object into mageck input.

        If screen.guides[target_column] is None or np.isnan, remove that guide
        """
        if count_layer is None:
            count_matrix = self.X
        else:
            try:
                count_matrix = self.layers[count_layer]
            except KeyError as exc:
                raise KeyError(
                    f"Layer {count_layer} doesn't exist in Screen object with layers {self.layers.keys()}"
                ) from exc
        mageck_input_df = (
            pd.DataFrame(
                count_matrix,
                columns=sample_prefix + self.condit.index,
                index=self.guides.index,
            )
            .fillna(0)
            .astype(int)
        )
        if sgrna_column is None:
            mageck_input_df.insert(0, "sgRNA", self.guides.index.tolist())
        elif sgrna_column in self.guides.columns:
            mageck_input_df.insert(0, "sgRNA", self.guides[sgrna_column])
        elif self.guides.index.name == sgrna_column:
            mageck_input_df.insert(0, "sgRNA", self.guides.index.tolist())
        else:
            raise ValueError(f"{sgrna_column} not found in Screen.guides.")
        mageck_input_df["sgRNA"] = mageck_input_df["sgRNA"].map(
            lambda s: s.replace(" ", "_")
        )
        mageck_input_df.insert(1, "gene", self.guides[target_column])
        mageck_input_df = mageck_input_df.loc[
            (mageck_input_df.gene.map(lambda o: not pd.isnull(o)))
            & mageck_input_df.gene.map(bool),
            :,
        ]
        if out_path is None:
            return mageck_input_df
        else:
            mageck_input_df.to_csv(out_path, sep="\t", index=False)

    def write(self, out_path):
        """
        Write .h5ad
        """
        self.obs = self.guides
        self.var = self.condit
        super().write(out_path)


def read_h5ad(filename):
    return _Screen.read_h5ad(filename)


def concat(screens, *args, **kwargs):
    adata = ad.concat(screens, *args, **kwargs)

    return _Screen(adata)


def read_csv(X_path=None, guide_path=None, condit_path=None, sep=","):
    if X_path is not None:
        X_df = pd.read_csv(X_path, delimiter=sep, header=0, index_col=0)
        X = X_df.values
    else:
        X = None
    guide_df = pd.read_csv(guide_path, sep=sep) if guide_path is not None else None
    condit_df = None if condit_path is None else pd.read_csv(condit_path, sep=sep)
    return _Screen(X=X, guides=guide_df, condit=condit_df)
