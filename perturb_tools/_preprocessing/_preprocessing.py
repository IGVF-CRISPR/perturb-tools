__module_name__ = "_ScreenModule.py"
__author__ = ", ".join(["Michael E. Vinyard", "Jayoung Kim Ryu"])
__email__ = ", ".join(["vinyard@g.harvard.edu", "jayoung_ryu@g.harvard.edu"])

import warnings
from typing import Optional, Union
import numpy as np
import pandas as pd
from anndata import AnnData
from .._arithmetic._funcs._log_fold_change import _log_fold_change
from ._supporting_functions._read_count_norm import _log_normalize_read_count


def log_norm(adata: AnnData, output_layer="lognorm_counts", read_count_layer=None):
    if read_count_layer is None:
        adata.layers[output_layer] = _log_normalize_read_count(adata.X)
    else:
        output_layer = f"lognorm_{read_count_layer}"
        adata.layers[output_layer] = _log_normalize_read_count(
            adata.layers[read_count_layer]
        )


def log_fold_change(
    adata: AnnData,
    sample1: str,
    sample2: str,
    lognorm_counts_key: str = "lognorm_counts",
    name: Optional[str] = None,
    out_guides_suffix: str = "lfc",
    return_result: bool = False,
):
    """
    Calculate log fold-change (LFC) between two samples across experimental conditions.

    This function calculates the LFC (sample1/sample2) using log-normalized values across experimental
    conditions for each gene in the dataset.

    Parameters
    ----------
    adata : `anndata.AnnData`
        Annotated data matrix containing the gene expression data.
    sample1 : `str`
        Name of the first sample to compare.
    sample2 : `str`
        Name of the second sample to compare.
    lognorm_counts_key : `str`, optional (default: "lognorm_counts")
        Key to access normalized count data in `adata`.
    name : `str`, optional (default: None)
        Name for the newly calculated LFC value. If not provided, the default is generated based on the input samples names.
    out_guides_suffix : `str`, optional (default: "lfc")
        Suffix string to append at the end of each guide sequence in the resulting calculation.
    return_result : `bool`, optional (default: False)
        Whether to return the result instead of storing it in `adata`.

    Returns
    -------
    lfc : `pandas.Series` or None
        Series containing LFC for each gene, or None if `return_result` is False.

    Raises
    ------
    ValueError
        If either of the specified samples are not present in the `adata.obs` DataFrame.
    """

    if "lognorm" not in lognorm_counts_key:
        warnings.warn(
            "The layer specified must be log-normalized values using screen.log_norm()."
        )

    if lognorm_counts_key not in adata.layers.keys():
        raise ValueError(
            "Specified normalized count isn't in your layer. First run screen.log_norm()."
        )

    sample1_idx = np.where(sample1 == adata.obs.index)[0]
    sample2_idx = np.where(sample2 == adata.obs.index)[0]
    if len(sample1_idx) != 1 or len(sample2_idx) != 1:
        if len(sample1_idx) == 0:
            print(f"No sampleition named {sample1} in Screen object.")
        else:
            print(f"Duplicate sampleition name {sample1} in Screen object")
        if len(sample2_idx) == 0:
            print(f"No sampleition named {sample2} in Screen object.")
        else:
            print(f"Duplicate sampleition name {sample2} in Screen object")
        raise ValueError("")

    lfc = _log_fold_change(
        adata.layers[lognorm_counts_key], sample1_idx[0], sample2_idx[0]
    )
    if return_result:
        return lfc
    else:
        adata.var[f"{sample1}_{sample2}.{out_guides_suffix}"] = lfc


def log_fold_change_reps(
    adata: AnnData,
    cond1: str,
    cond2: str,
    lognorm_counts_key: str = "lognorm_counts",
    rep_col: str = "replicate",
    compare_col: str = "sort",
    out_guides_suffix: str = "lfc",
    keep_result: bool = False,
):
    if rep_col not in adata.obs.columns:
        raise ValueError(f"{rep_col} not in condit features")
    if compare_col not in adata.obs.columns:
        raise ValueError(f"{compare_col} not in condit features")

    lfcs = []
    for rep in adata.obs[rep_col].unique():
        cond1_idx = np.where(
            (adata.obs[rep_col] == rep) & (adata.obs[compare_col] == cond1)
        )[0]
        cond2_idx = np.where(
            (adata.obs[rep_col] == rep) & (adata.obs[compare_col] == cond2)
        )[0]

        if len(cond1_idx) != 1 or len(cond2_idx) != 1:
            raise ValueError(
                "Conditions are not unique for each replicates to be aggregated."
            )

        lfcs.append(
            log_fold_change(
                adata,
                adata.obs.index[cond1_idx].tolist()[0],
                adata.obs.index[cond2_idx].tolist()[0],
                lognorm_counts_key=lognorm_counts_key,
                return_result=True,
            )
        )

    lfcs_array = np.vstack(lfcs).T
    lfcs_df_columns = [
        f"{s}.{cond1}_{cond2}.{out_guides_suffix}" for s in adata.obs[rep_col].unique()
    ]
    lfcs_df = pd.DataFrame(lfcs_array, index=adata.var.index, columns=lfcs_df_columns)

    if keep_result:
        adata.var[lfcs_df_columns] = lfcs_df
    return lfcs_df


# TODO: add guides metadata on how aggregates are calcualted?
def log_fold_change_aggregate(
    adata: AnnData,
    cond1: str,
    cond2: str,
    lognorm_counts_key: str = "lognorm_counts",
    aggregate_col: str = "replicate",
    compare_col: str = "sort",
    out_guides_suffix: str = "lfc",
    aggregate_fn: str = "median",
    name: Optional[str] = None,
    return_result: bool = False,
    keep_per_replicate: bool = False,
):
    """Calculates aggregated log fold changes using specified aggregation method and stores the result as metadata in AnnData object.

    Parameters:
    ----------
    adata: `anndata.AnnData`
        Annotated data matrix containing the gene expression data.
    cond1: `str`
        Name of the first condition to compare.
    cond2 : `str`
        Name of the second condition to compare.
    lognorm_counts_key: `str`, optional (default: "lognorm_counts")
        Key for accessing normalized count data in `adata`.
    aggregate_col: `str`, optional (default: "replicate")
        Names of conditions to aggregate over.
    compare_col: `str`, optional (default: "sort")
        Names of conditions to perform comparison between.
    out_guides_suffix: `str`, optional (default: "lfc")
        Suffix for guide.
    aggregate_fn: `str`, optional (default: "median")
        Method for aggregating log fold changes. Supported methods are "mean", "median" and "sd".
    name: `str`, optional (default: None)
        Name of metadata column to store result.
    return_result: `bool`, optional (default: False)
        Whether to return the result instead of storing it in `adata`.
    keep_per_replicate: `bool`, optional (default: False)
        Whether to keep individual replicate log fold changes.

    Returns:
    --------
    lfcs_agg: `pandas.DataFrame` or None
        Dataframe containing aggregated log fold changes, or None if `return_result` is False.

    Raises:
    -------
    ValueError
        If the provided `aggregate_fn` parameter is not "mean", "median", or "sd".

    """
    lfcs_df = log_fold_change_reps(
        adata,
        cond1,
        cond2,
        lognorm_counts_key=lognorm_counts_key,
        rep_col=aggregate_col,
        compare_col=compare_col,
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
        adata.var[f"{cond1}_{cond2}.{out_guides_suffix}.{aggregate_fn}"] = lfcs_agg
    else:
        adata.var[name] = lfcs_agg


def fold_change(
    adata: AnnData,
    cond1: str,
    cond2: str,
    lognorm_counts_key: str = "lognorm_counts",
    return_result: bool = False,
) -> Union[None, pd.Series]:
    """Calculate the log-fold change between two conditions.

    Parameters
    ----------
    adata : `anndata.AnnData`
        Annotated data matrix containing the gene expression data.
    cond1 : `str`
        Name of the first condition to compare.
    cond2 : `str`
        Name of the second condition to compare.
    lognorm_counts_key : `str`, optional (default: "lognorm_counts")
        Key to access normalized count data in `adata`.
    return_result : `bool`, optional (default: False)
        Whether to return the result instead of storing it in `adata`.

    Returns
    -------
    lfcs : `pandas.Series` or None
        Series containing log-fold changes for each gene, or None if `return_result` is False.

    Raises
    ------
    ValueError
        If either of the specified conditions are not present in the `adata.obs` DataFrame.
    """

    log_fold_change = _log_fold_change(adata.layers[lognorm_counts_key], cond1, cond2)
    if return_result:
        return log_fold_change
    adata.var[f"{cond1}_{cond2}.fc"] = log_fold_change
