import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def _calculate_enrichment_pvalues(arr1, arr2):

    """
    Calculate the probability that two values drawn N times (len(row)) are from distinct distributions.


    Calculate the T-test for the means of two independent samples of scores.
    (source: scipy.stats)


    Notes:
    ------
    (1) This is a two-sided test for the null hypothesis that 2
        independent samples have identical average (expected)
        values.

    (2) This test assumes that the populations have identical
        variances by default.

    (3) Notes (1) and (2) are taken directly from the scipy.stats website.
    """

    p_vals = []

    for row in range(len(arr1)):

        cond1 = arr1[row]
        cond2 = arr2[row]

        p_vals.append(scipy.stats.ttest_ind(cond1, cond2)[1])

    return np.array(p_vals)


def _calculate_baseline_subtraction(condit_1, condit_2, control):

    """Calculate the log fold change relative to the control (i.e., perform baseline subtraction).

    Assumes log-transformed values."""

    condit_1_vs_control = condit_1.values - control.values
    condit_2_vs_control = condit_2.values - control.values

    return condit_1_vs_control, condit_2_vs_control


def _calculate_delta_logfoldchange(condit_1, condit_2):

    """
    Calculate the delta log fold change in guide counts between
    two conditions. Calculates the mean and standard deviation
    over replicates.

    Parameters:
    -----------
    condit_1

    condit_2

    Returns:
    --------


    Notes:
    ------
    (1) Assumes baseline has already been substracted.

    """

    log_fold_change_mean = condit_1.mean(axis=1) - condit_2.mean(axis=1)
    log_fold_change_stdev = (condit_1 - condit_2).std(axis=1)

    return log_fold_change_mean, log_fold_change_stdev