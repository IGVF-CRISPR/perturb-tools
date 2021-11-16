
# _ScreenModule.py
__module_name__ = "_ScreenModule.py"
__author__ = ", ".join(["Michael E. Vinyard", "Jayoung Kim Ryu"])
__email__ = ", ".join(["vinyard@g.harvard.edu", "jayoung_ryu@g.harvard.edu"])


import pandas as pd

from ._supporting_functions._print_screen_object import _print_screen_object
from ._supporting_functions._data_reading._read_screen_from_PoolQ import _read_screen_from_PoolQ
from ._supporting_functions._guides._GuideAnnotationModule import _annotate_sgRNAs
from .._normalization._funcs._read_count_norm import _log_normalize_read_count
from ._supporting_functions._print_screen_object import _print_screen_object

from .._arithmetic._funcs._log_fold_change import _log_fold_change
from .._arithmetic._funcs._fold_change import _fold_change

from .._readwrite._funcs._write_screen_to_csv import _write_screen_to_csv
from .._readwrite._funcs._write_screen_to_excel import _write_screen_to_excel

from .._utilities._funcs._update_dict import _update_dict

class _Screen:
    def __init__(self, X=None, *args, **kwargs):
        if X is not None:
            ad = AnnData(X, *args, **kwargs)
            self.X = ad.X
            self.guides = ad.obs
            self.condit = ad.var
            self.condit_m = ad.obsm
            self.condit_p = ad.obsp
            self.layers = ad.layers
            self.uns = ad.uns
            n_guides, n_conditions, _ = _print_screen_object(self)
            
    def __repr__(self) -> str:
        return _print_screen_object(self)[2]

    def read_PoolQ(self, path, metadata=False, merge_metadata_on='Condition'):

        """ 
        Read poolQ.
        """

        self._PoolQ_outpath = path
        self._PoolQScreenDict = _read_screen_from_PoolQ(self._PoolQ_outpath)
        
        for key, value in v.ut.update_dict(self._PoolQScreenDict).items():
            self.__setattr__(key, value)
        
        if metadata:
            self.condit = self.condit.merge(pd.read_csv(metadata), on=merge_metadata_on)
            
        _print_screen_object(self)
        
        
    def annotate_guides(self, genes, chrom, start, stop, annotations, DirectPairDict, ref_seq_path):
        
        """
        Annotate sgRNA table.

        """
        self.guides = _annotate_sgRNAs(self.guides, genes, chrom, start, stop, annotations, DirectPairDict, ref_seq_path)
        
    def log_norm(self, layer_key='lognorm_counts'):
        
        self.layers[layer_key] = _log_normalize_read_count(self.X)
        
    
    def log_fold_change(self, cond1, cond2, lognorm_counts_key="lognorm_counts", name=False):

        """ 
        General module to calculate LFC across experimental conditions. 
        """
        try:
            self.guides["{}_{}.lfc".format(cond1, cond2)] = _log_fold_change(
                self.layers[lognorm_counts_key], cond1, cond2
            )
        except:
            print("Calculating LFC against two previously calculated LFC values...")
            
            if not name:
                name = "{}_{}.dlfc".format(cond1.strip(".lfc"), cond2.strip(".lfc"))
            
            self.guides[name] = _log_fold_change(
                self.guides, cond1, cond2
            )
        
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
