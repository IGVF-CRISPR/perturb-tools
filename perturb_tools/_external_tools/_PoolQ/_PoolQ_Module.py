# package imports #
# --------------- #
import numpy as np
import os

# local imports #
# ------------- #
from ._supporting_functions._PoolQ_supporting_functions import _get_read_files
from ._supporting_functions._PoolQ_supporting_functions import _run_PoolQ
from ._supporting_functions._PoolQ_supporting_functions import _organize_PoolQ_outputs

# from ._PoolQ_supporting_functions import _return_outs_filepaths
from ._supporting_functions._PoolQ_supporting_functions import _make_results_dict
from ._supporting_functions._PoolQ_supporting_functions import _plot_correlation_scatter
from ._supporting_functions._get_guide_position_df import _get_guide_position_df
from ._supporting_functions._logfoldchange import _plot_multiLFC

from ..._utilities._funcs._flexible_mkdir import _flexible_multilevel_mkdir
from ..._plotting._funcs._plot_correlation_heatmap import _plot_correlation_heatmap

poolq_path = os.path.join(__file__, "_external_software/poolq-3.3.2/poolq3.sh")

class _PoolQ:
    def __init__(
        self,
        sequencing,
        metadata,
        outpath,
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
        data_dir,
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
            data_dir,
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

    def plot_correlation(
        self,
        title="Sample Correlation",
        figsize=3,
        title_y=1.15,
        title_x=-0.1,
        save=False,
    ):

        """
        
        .png already included in savename. 
        """

        if save:
            savename = os.path.join(
                self.outpath, self.run_name + ".correlation.heatmap"
            )
        else:
            savename = False

        figsavename = _plot_correlation_heatmap(
            self.ResultsDict["correlation"],
            title=title,
            figsize=figsize,
            title_y=title_y,
            title_x=title_x,
            savename=savename,
        )

        print("\nHeatmap saved to: {}".format(figsavename))

        if save:
            savename = os.path.join(self.outpath, self.run_name)
        else:
            savename = False

        self.corr_fig, figsavename = _plot_correlation_scatter(
            self.ResultsDict["lognorm_counts"],
            savename=savename,
            figsize=figsize,
            ignore_containing=["sequence", "barcode_id"],
        )

        print("\nCorrelation scatter plot saved to: {}".format(figsavename))

    def get_guide_positions(self, TargetDict, row_delim="\t"):

        region = TargetDict["region"]
        chromosome = TargetDict["chromosome"]
        start, stop = TargetDict["start"], TargetDict["stop"]
        self.rowfile = TargetDict["rowsfile"]

        self.guide_df = _get_guide_position_df(
            region, chromosome, start, stop, self.rowfile, row_delim
        )

    def log_fold_change(self, conditions, controls, ms=25, elinewidth=4, narrow=True):

        _plot_multiLFC(
            df=self.ResultsDict["lognorm_counts"],
            rowfile=self.rowfile,
            guide_df=self.guide_df,
            conditions=conditions,
            controls=controls,
            ms=ms,
            elinewidth=elinewidth,
            narrow=narrow,
        )
