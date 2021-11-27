# package imports #
# --------------- #
import numpy as np, pandas as pd
import os
import string

# local imports #
# ------------- #
from ...._utilities._funcs._flexible_mkdir import _flexible_multilevel_mkdir
from ._supporting_functions._run_PoolQ import _run_PoolQ


def _read_txt(filepath):

    """
    Read a text file given a filepath.
    
    Parameters:
    -----------
    filepath
        path to textfile
        type:str
    
    Returns:
    --------
    lines
        lines of content in the file appended into a list
        type: list
    
    Notes:
    ------
    (1) Skips lines without content. 
    """

    lines = []

    f = open(filepath, "r")
    for line in f.readlines():
        if len(line) == 1:
            continue
        else:
            lines.append(line.strip("\n"))

    return lines


def _make_results_dict(poolq_outs_path):

    ResultsDict = {}

    ResultsDict["raw_counts"] = pd.read_csv(
        os.path.join(poolq_outs_path, "counts.txt"), sep="\t"
    )
    ResultsDict["barcode_counts"] = pd.read_csv(
        os.path.join(poolq_outs_path, "barcode-counts.txt"), sep="\t"
    )
    ResultsDict["correlation"] = pd.read_csv(
        os.path.join(poolq_outs_path, "correlation.txt"), sep="\t", index_col=0
    )
    ResultsDict["lognorm_counts"] = pd.read_csv(
        os.path.join(poolq_outs_path, "lognormalized-counts.txt"), sep="\t"
    )

    ResultsDict["poolq3.log"] = _read_txt(os.path.join(poolq_outs_path, "poolq3.log"))
    ResultsDict["quality"] = _read_txt(os.path.join(poolq_outs_path, "quality.txt"))
    ResultsDict["runinfo"] = _read_txt(os.path.join(poolq_outs_path, "runinfo.txt"))
    ResultsDict["unexpected-sequences"] = pd.read_csv(
        os.path.join(poolq_outs_path, "unexpected-sequences.txt"), sep="\t"
    )

    return ResultsDict



def _organize_PoolQ_outputs(run_name, standard_PoolQ_outfiles, outpath, fetch_only):

    """Move the outputs of PoolQ from the default output directory (i.e., the working directory "./") to one defined."""

    run_dir = run_name + ".PoolQ.3.3.2"

    if not fetch_only:
        _flexible_multilevel_mkdir(run_dir)

    for file in standard_PoolQ_outfiles:
        if not fetch_only:
            os.system("mv {} {}".format(file, os.path.join(run_dir)))

    if outpath:
        if not fetch_only:
            os.system("mv {} {}".format(run_dir, outpath))
        outs = os.path.join(outpath, run_dir)
        return outs
    else:
        outs = os.path.join(os.getcwd(), run_dir)
        return outs



# from ._supporting_functions._PoolQ_supporting_functions import _get_read_files
# from ._supporting_functions._PoolQ_supporting_functions import _run_PoolQ
# from ._supporting_functions._PoolQ_supporting_functions import _organize_PoolQ_outputs

# # from ._PoolQ_supporting_functions import _return_outs_filepaths
# from ._supporting_functions._PoolQ_supporting_functions import _make_results_dict
# from ._supporting_functions._PoolQ_supporting_functions import _plot_correlation_scatter
# from ._supporting_functions._get_guide_position_df import _get_guide_position_df
# from ._supporting_functions._logfoldchange import _plot_multiLFC

# from ...._utilities._funcs._flexible_mkdir import _flexible_multilevel_mkdir
# from ...._plotting._funcs._plot_correlation_heatmap import _plot_correlation_heatmap

poolq_path = os.path.join(os.path.dirname(__file__), "_poolq_software/poolq-3.3.2/poolq3.sh")

class _PoolQ:
    def __init__(
        self,
        data_dir=False,
        metadata=False,
        outpath=False,
        verbosity=True,
        barcode_indicator="barcode",
        barcode_filename="rows.txt",
        conditions_filename="column.txt",
        standard_PoolQ_outfiles=[
            "barcode-counts.txt",
            "correlation.txt",
            "counts.txt",
            "lognormalized-counts.txt",
            "poolq3.log",
            "quality.txt",
            "runinfo.txt",
            "unexpected-sequences.txt",
        ],
    ):

        self.data_dir = data_dir
        self.barcode_indicator = barcode_indicator
        self.barcode_filename = barcode_filename
        self.conditions_filename = conditions_filename
        self.standard_PoolQ_outfiles = standard_PoolQ_outfiles
        if outpath:
            _flexible_multilevel_mkdir(outpath, verbose=verbosity)
            self.outpath = outpath

    #             self.outfiles = _return_outs_filepaths(self.outpath, self.standard_PoolQ_outfiles)

    def run(
        self,
        data_dir=False,
        run_name=False,
        dry_run=False,
        barcode_filename=False,
        conditions_filename=False,
        barcode_indicator=False,
        outpath=False,
    ):
        if not run_name:
            self.run_name = "run_{}".format(
                "".join(np.random.choice(list(string.ascii_lowercase), 6))
            )
        else:
            self.run_name = run_name

        if data_dir:
            self.data_dir = data_dir
        if barcode_filename:
            self.barcode_filename = barcode_filename
        if conditions_filename:
            self.conditions_filename = conditions_filename
        if barcode_indicator:
            self.barcode_indicator = barcode_indicator
        if outpath:
            self.outpath = outpath

        _run_PoolQ(
            self.data_dir,
            self.barcode_filename,
            self.conditions_filename,
            poolq_path,
            self.barcode_indicator,
            dry_run,
            self.run_name,
        )

        if not dry_run:
            self.outs = _organize_PoolQ_outputs(
                self.run_name,
                self.standard_PoolQ_outfiles,
                self.outpath,
                fetch_only=False,
            )
            self.ResultsDict = _make_results_dict(self.outs)

    def load_results(
        self,
        run_name=False,
        standard_PoolQ_outfiles=False,
        outpath=False,
        fetch_only=False,
    ):

        """"""

        if run_name:
            self.run_name = run_name

        if standard_PoolQ_outfiles:
            self.standard_PoolQ_outfiles = standard_PoolQ_outfiles

        if outpath:
            self.outpath = outpath

        self.outs = _organize_PoolQ_outputs(
            self.run_name, self.standard_PoolQ_outfiles, self.outpath, fetch_only
        )
        self.ResultsDict = _make_results_dict(self.outs)

#     def plot_correlation(
#         self,
#         title="Sample Correlation",
#         figsize=3,
#         title_y=1.15,
#         title_x=-0.1,
#         save=False,
#     ):

#         """
        
#         .png already included in savename. 
#         """

#         if save:
#             savename = os.path.join(
#                 self.outpath, self.run_name + ".correlation.heatmap"
#             )
#         else:
#             savename = False

#         figsavename = _plot_correlation_heatmap(
#             self.ResultsDict["correlation"],
#             title=title,
#             figsize=figsize,
#             title_y=title_y,
#             title_x=title_x,
#             savename=savename,
#         )

#         print("\nHeatmap saved to: {}".format(figsavename))

#         if save:
#             savename = os.path.join(self.outpath, self.run_name)
#         else:
#             savename = False

#         self.corr_fig, figsavename = _plot_correlation_scatter(
#             self.ResultsDict["lognorm_counts"],
#             savename=savename,
#             figsize=figsize,
#             ignore_containing=["sequence", "barcode_id"],
#         )

#         print("\nCorrelation scatter plot saved to: {}".format(figsavename))

#     def get_guide_positions(self, TargetDict, row_delim="\t"):

#         region = TargetDict["region"]
#         chromosome = TargetDict["chromosome"]
#         start, stop = TargetDict["start"], TargetDict["stop"]
#         self.rowfile = TargetDict["rowsfile"]

#         self.guide_df = _get_guide_position_df(
#             region, chromosome, start, stop, self.rowfile, row_delim
#         )

#     def log_fold_change(self, conditions, controls, ms=25, elinewidth=4, narrow=True):

#         _plot_multiLFC(
#             df=self.ResultsDict["lognorm_counts"],
#             rowfile=self.rowfile,
#             guide_df=self.guide_df,
#             conditions=conditions,
#             controls=controls,
#             ms=ms,
#             elinewidth=elinewidth,
#             narrow=narrow,
#         )
