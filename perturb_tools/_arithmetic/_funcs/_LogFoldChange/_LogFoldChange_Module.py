
import numpy as np

from ._calculate_log_fold_change import _calculate_enrichment_pvalues
from ._calculate_log_fold_change import _calculate_baseline_subtraction
from ._calculate_log_fold_change import _calculate_delta_logfoldchange

#from ...._plotting._funcs._plot_guide_enrichment import _plot_guide_enrichment

class LogFoldChangeModule:
    def __init__(self, counts_df):

        """
        Calculate the log fold change between two conditions.

        Parameters:
        -----------
        counts_df
            example: screen.layers["lognorm_counts"]
            type: pandas.core.frame.DataFrame

        condit_1
            Regex string by which to identify condition columns in counts_df that match condit_1.
            type: str

        condit_2
            Regex string by which to identify condition columns in counts_df that match condit_2.
            type: str

        control
            Regex string by which to identify condition columns in counts_df that match control condition.
            type: str


        Returns:
        --------

        Notes:
        ------
        (1) Assumes that two conditions are to be compared, originating from a common initial condition.
        (2) Assumes log-transformed values.
        """

        self.counts_df = counts_df

    def filter_guides_by_experiments(self, guides, targets):

        """
        Remove guides that do not target the loci of interest.

        Control guides may be shared across experiments.

        Parameters:
        -----------
        counts_df

        screen.guides

        list(str)

        Notes:
        ------
        counts_df and guides should have the same length with (initially) matching indices.
        """

        self.counts_df = self.counts_df.iloc[
            guides["target"].loc[guides["target"].isin(targets)].index.astype(int)
        ]

    def isolate_conditions(self, condit_1, condit_2, control):

        self.condit_1_ = self.counts_df.filter(regex=condit_1)
        self.condit_2_ = self.counts_df.filter(regex=condit_2)
        self.control_ = self.counts_df.filter(regex=control)

    def calculate_logfoldchange(self):

        """
        Three steps:
        ------------
        (1) calculate baseline.
        (2) calculate baseline-LFC.
        (3) Calculate a p-value for the logfoldchange of each guide between conditions.

        """

        (
            self.condit_1_vs_control,
            self.condit_2_vs_control,
        ) = _calculate_baseline_subtraction(
            self.condit_1_, self.condit_2_, self.control_
        )

        self.lfc_mean, self.lfc_stdev = _calculate_delta_logfoldchange(
            self.condit_1_vs_control, self.condit_2_vs_control
        )

        self.p_values = _calculate_enrichment_pvalues(
            self.condit_1_vs_control, self.condit_2_vs_control
        )
        
        
def _annotate_guides_with_log_fold_change(screen, lfc_module):

    """"""

    experiment_guide_df = screen.guides.iloc[lfc_module.counts_df.index]
    experiment_guide_df["lfc.mean"] = lfc_module.lfc_mean
    try:
        experiment_guide_df["lfc.stdev"] = lfc_module.lfc_stdev
        experiment_guide_df["lfc.p_value"] = lfc_module.p_values
        experiment_guide_df["-log10(pval)"] = -np.log10(experiment_guide_df["lfc.p_value"])
    except:
        pass

    return experiment_guide_df
        
    
def _log_fold_change_analysis(
    screen,
    condit_1,
    condit_2,
    control,
    targets,
    layer="lognorm_counts",
    plot=True,
    plot_title=None,
    return_fig=False,
    gene=False,
    gtf_path=None,
):

    """"""

    lfc = LogFoldChangeModule(screen.layers[layer])
    lfc.filter_guides_by_experiments(screen.guides, targets)
    lfc.isolate_conditions(condit_1, condit_2, control)
    lfc.calculate_logfoldchange()
    
    experiment_guide_df = _annotate_guides_with_log_fold_change(screen, lfc)
    
    if plot:
        fig = _plot_guide_enrichment(experiment_guide_df, title=plot_title, return_fig=return_fig, gene=gene, gtf_path=gtf_path)
        if return_fig:
            return experiment_guide_df, fig
    else:
        return experiment_guide_df
