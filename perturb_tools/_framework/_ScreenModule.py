
# _ScreenModule.py
__module_name__ = "_ScreenModule.py"
__author__ = ", ".join(["Michael E. Vinyard", "Jayoung Kim Ryu"])
__email__ = ", ".join(["vinyard@g.harvard.edu", "jayoung_ryu@g.harvard.edu"])


import pandas as pd
import numpy as np
import anndata as ad
import copy
import warnings

from ._supporting_functions._print_screen_object import _print_screen_object
from ._supporting_functions._guides._GuideAnnotationModule import _annotate_sgRNAs
from .._normalization._funcs._read_count_norm import _log_normalize_read_count
from ._supporting_functions._print_screen_object import _print_screen_object

from .._readwrite._funcs._write_screen_to_csv import _write_screen_to_csv
from .._readwrite._funcs._write_screen_to_excel import _write_screen_to_excel
from .._readwrite._funcs._read_screen_from_PoolQ import _read_screen_from_PoolQ

from .._utilities._funcs._update_dict import _update_dict

class _Screen:
    ad = __import__('anndata')
    def __init__(self, X=None, *args, **kwargs):
        if X is not None:
            adata = ad.AnnData(X, *args, **kwargs)
            self.X = adata.X
            self.guides = adata.obs
            self.condit = adata.var
            self.condit_m = adata.obsm
            self.condit_p = adata.obsp
            self.layers = adata.layers
            self.uns = adata.uns
            n_guides, n_conditions, _ = _print_screen_object(self)
            
    def __repr__(self) -> str:
        return _print_screen_object(self)[2]


    def __add__(self, other):
        if all(self.guides.index == other.guides.index) and all(self.condit.index == other.condit.index):
            added = self.copy()
            for layer in self.layers:
                added.layers[layer] += other.layers[layer]
            return(added)
        else:
            raise ValueError("Guides/sample description mismatch")

    def copy(self):
        return copy.deepcopy(self)


    def read_PoolQ(self, path, metadata=False, merge_metadata_on='Condition'):

        """ 
        Read poolQ.
        """
        
        self._PoolQ_outpath = path
        self._PoolQScreenDict = _read_screen_from_PoolQ(self._PoolQ_outpath)
        
        for key, value in _update_dict(self._PoolQScreenDict).items():
            self.__setattr__(key, value)
        
        if metadata:
            self.condit = self.condit.merge(pd.read_csv(metadata), on=merge_metadata_on)
            
        _print_screen_object(self)
        
        
    def annotate_guides(self, genes, chrom, start, stop, annotations, DirectPairDict, ref_seq_path):
        
        """
        Annotate sgRNA table.

        """
        self.guides = _annotate_sgRNAs(self.guides, genes, chrom, start, stop, annotations, DirectPairDict, ref_seq_path)
        
    def log_norm(self, output_layer='lognorm_counts', read_count_layer = None):
        if read_count_layer is None:
            self.layers[output_layer] = _log_normalize_read_count(self.X)
        else:
            output_layer = "lognorm_{}".format(read_count_layer)
            self.layers[output_layer] = _log_normalize_read_count(self.layers[read_count_layer])
        
    
    # TBD: mask ones with too low raw counts.
    def log_fold_change(
            self, 
            cond1, 
            cond2, 
            lognorm_counts_key="lognorm_counts", 
            name=False, 
            out_guides_suffix = "lfc", 
            return_result = False
            ):

        """ 
        General module to calculate LFC across experimental conditions. 
        """
        if "lognorm" not in lognorm_counts_key:
            warnings.warn("The layer specified must be log-normalized values using screen.log_norm().")

        if lognorm_counts_key not in self.layers.keys(): 
            raise ValueError("Specified normalized count isn't in your layer. First run screen.log_norm().")

        
        cond1_idx = np.where(cond1 == self.condit.index)[0]
        cond2_idx = np.where(cond2 == self.condit.index)[0]
        if len(cond1_idx) != 1 or len(cond2_idx) != 1:
            if len(cond1_idx) == 0: print("No condition named {} in Screen object.".format(cond1))
            else: print("Duplicate condition name {} in Screen object".format(cond1))
            if len(cond2_idx) == 0: print("No condition named {} in Screen object.".format(cond2))
            else: print("Duplicate condition name {} in Screen object".format(cond2))
            raise ValueError("")


        try:
            lfc = _log_fold_change(
                self.layers[lognorm_counts_key], cond1_idx, cond2_idx
            )
            if return_result:
                return(lfc)
            else:
                self.guides["{}_{}.{}".format(cond1, cond2, out_guides_suffix)] = lfc
        except:# TBD: what error?
            print("Calculating LFC against two previously calculated LFC values...")
            dlfc = _log_fold_change(self.guides, cond1, cond2)
            
            if not name:
                name = "{}_{}.d{}".format(cond1.strip(".lfc"), cond2.strip(".lfc"), out_guides_suffix)
            
            if return_result:
                return(dlfc)
            else:
                self.guides[name] = dlfc

    def log_fold_change_reps(
            self, 
            cond1, 
            cond2,
            lognorm_counts_key = "lognorm_counts",
            aggregate_condit = "replicate", 
            compare_condit = "sort",
            out_guides_suffix = "lfc",
            keep_result = False,
            ):
        
        if aggregate_condit not in self.condit.columns:
            raise ValueError("{} not in condit features".format(aggreagate_condit))
        if compare_condit not in self.condit.columns:
            raise ValueError("{} not in condit features".format(compare_condit))


        lfcs = []
        for rep in self.condit[aggregate_condit].unique():
            cond1_idx = np.where((self.condit[aggregate_condit] == rep) & (self.condit[compare_condit] == cond1))[0]
            cond2_idx = np.where((self.condit[aggregate_condit] == rep) & (self.condit[compare_condit] == cond2))[0]

            if len(cond1_idx) != 1 or len(cond2_idx) != 1:
                raise ValueError("Conditions are not unique for each replicates to be aggregated.")
            
            lfcs.append(self.log_fold_change(
                self.condit.index[cond1_idx].tolist()[0], 
                self.condit.index[cond2_idx].tolist()[0],
                lognorm_counts_key = lognorm_counts_key,
                return_result = True
                ))

        lfcs_array = np.concatenate(lfcs, axis = 1)
        lfcs_df_columns =  [s + ".{}_{}.{}".format(cond1, cond2, out_guides_suffix) for s in self.condit[aggregate_condit].unique()] 
        lfcs_df = pd.DataFrame(lfcs_array, 
                index = self.guides.index,
                columns = lfcs_df_columns )

        if keep_result : self.guides[lfcs_df_columns] = lfcs_df
        return(lfcs_df)
                        

    # TODO: add guides metadata on how aggregates are calcualted?
    def log_fold_change_aggregate(
            self, 
            cond1, 
            cond2,
            lognorm_counts_key = "lognorm_counts",
            aggregate_condit = "replicate", 
            compare_condit = "sort",
            out_guides_suffix = "lfc",
            aggregate_fn = "median",
            name = None,
            return_result = False,
            keep_per_replicate = False,
            ):
        
        lfcs_df = self.log_fold_change_reps(
                cond1, 
                cond2,
                lognorm_counts_key = lognorm_counts_key,
                aggregate_condit = aggregate_condit, 
                compare_condit = compare_condit,
                out_guides_suffix = out_guides_suffix,
                keep_result = keep_per_replicate,
                )

        if aggregate_fn == "mean":
            lfcs_agg = lfcs_df.apply(np.mean, axis = 1)
        elif aggregate_fn == "median":
            lfcs_agg = lfcs_df.apply(np.median, axis = 1)
        elif aggregate_fn == "sd":
            lfcs_agg = lfcs_df.apply(np.std, axis = 1)
        else:
            raise ValueError("Only 'mean', 'median', and 'sd' are supported for aggregating LFCs.")

        if return_result: 
            return(lfcs_agg)
        else:
            if name is None:
                self.guides["{}_{}.{}.{}".format(cond1, cond2, out_guides_suffix, aggregate_fn)] = lfcs_agg
            else:
                self.guides[name] = lfcs_agg
            



    def fold_change(self, cond1, cond2):

        """
        # incomplete
        """

        self.guides["{}_{}.fc".format(cond1, cond2)] = _fold_change(
            self.layers[lognorm_counts_key], cond1, cond2
        )
        
    def to_Excel(self, workbook_path="CRISPR_screen.workbook.xlsx", index=False, silent=False):
        
        """
        Write components of Screen class to an Excel workbook. 
        
        Parameters:
        -----------
        workbook_path
            Prevent printing outpaths / details of the created workbook. 
            default: "CRISPR_screen.workbook.xlsx"
            type: str
            
        index
            Include an index in the workbook sheet. 
            default: False
            type: bool
            
        silent
            Prevent printing outpaths / details of the created workbook. 
            default: False
            type: bool
        
        Returns:
        --------
        None, writes to excel workbook.
        
        Notes:
        ------
        (1) Will likely need to be updated once we fully transition over to AnnData-like class. 
        
        """
                
        _write_screen_to_excel(self,
                                 workbook_path,
                                 index,
                                 silent,)
        
    def to_csv(self, out_path="CRISPR_screen"):
        
        """
        
        Write .csv files for each part of the screen. will eventually be replaced by something more native to AnnData. 
        
        """
        
        
        _write_screen(self, out_path)
