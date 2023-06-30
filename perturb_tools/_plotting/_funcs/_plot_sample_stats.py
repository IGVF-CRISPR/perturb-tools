from typing import Tuple, Optional, Union, List
import numpy as np
import pandas as pd
from anndata import AnnData
import matplotlib.axes
import matplotlib.pyplot as plt
import seaborn as sns
from ..._qc.qc import _set_sample_correlation
from ..._preprocessing._preprocessing import log_fold_change_reps


def _G(v: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Calculate the Gini coefficient for a given numpy array.

    Args:
        v (np.ndarray): Array for which the Gini coefficient is to be calculated.

    Returns:
        Tuple[np.ndarray, np.ndarray, float]: A tuple of three elements containing numpy arrays:
            1. An array of quantiles used for calculating Gini
            2. An array of corresponding bin fractions
            3. A float indicating the Gini coefficient

    References:
        - https://stackoverflow.com/questions/39512260/calculating-gini-coefficient-in-python-numpy
        - http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    """
    bins = np.linspace(0.0, 1.0, 50)
    total = float(np.sum(v))
    yvals = []
    for b in bins:
        bin_vals = v[v <= np.quantile(v, b)]
        bin_fraction = np.sum(bin_vals) / total
        yvals.append(bin_fraction)
    # perfect equality area
    pe_area = np.trapz(bins, x=bins)
    # lorenz area
    lorenz_area = np.trapz(yvals, x=bins)
    gini_val = (pe_area - lorenz_area) / float(pe_area)
    return bins, yvals, gini_val


def sample_count_gini(
    adata: AnnData,
    ax=None,
    figsize=(6, 4),
    count_layer: Optional[str] = None,
    *args,
    **kwargs,
):
    """Plot Gini index of guide coverage for each sample"""
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)
    if not count_layer:
        count_matrix = adata.X
    else:
        count_matrix = adata.layers[count_layer]
    for i, c in enumerate(adata.obs.index):
        bins, result, gini_val = _G(count_matrix[i, :])
        ax.plot(bins, result, label="{} ({:.3f})".format(c, gini_val), *args, **kwargs)
    ax.plot((0, 1), (0, 1), "--", label="perfect eq.")
    ax.set_xlabel("Cumulative fraction of guides \nfrom lowest to highest read counts")
    ax.set_ylabel("Cumulative fraction of read counts")
    ax.set_title("GINI index")
    ax.legend(bbox_to_anchor=(1.02, 1))


def sample_count_dist(
    adata: AnnData,
    ax=None,
    figsize=(6, 4),
    count_layer: Optional[str] = None,
    n_bins: int = 100,
    log_x: bool = True,
    *args,
    **kwargs,
):
    """Plot guide count distribution per sample"""
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)
    if count_layer is None:
        count_matrix = adata.X
    else:
        count_matrix = adata.layers[count_layer]
    if log_x:
        bins = 10 ** np.linspace(0, np.log10(count_matrix.max() + 1), n_bins)
    else:
        bins = np.linspace(0, count_matrix.max(), n_bins)
    for i, c in enumerate(adata.obs.index):
        if count_layer is None:
            count_matrix = adata.X[i, :]
        else:
            count_matrix = adata.layers[count_layer][i, :]
        ax.hist(
            count_matrix,
            bins,
            histtype="step",
            label=f"{c} (median={np.nanmedian(count_matrix):.3g})",
            *args,
            **kwargs,
        )
    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel("# guides")
    ax.set_ylabel("Count")
    ax.legend(bbox_to_anchor=(1.02, 1))
    return ax


def sample_count_correlation(
    screen: AnnData,
    count_layer: Optional[str] = None,
    method: str = "spearman",
    figsize: Tuple[int, int] = (10, 10),
    ax: Optional[matplotlib.axes.Axes] = None,
    guide_idx: Optional[Union[np.ndarray, List[bool]]] = None,
    correlation_prefix: str = "",
    *args,
    **kwargs,
) -> None:
    """Plot heatmap of pairwise correlation between samples in AnnData object based on guide counts.

    Args:
        screen (AnnData): Annotated data matrix containing guide counts.
        count_layer (str, optional): Key for the layer to use if AnnData object contains multiple layers. Defaults to None.
        method (str, optional): Method used for calculating correlation. Should be one of ["pearson", "spearman"].
        figsize ((int, int), optional): Width and height of figure in inches. Defaults to (10, 10).
        ax (matplotlib.axes.Axes, optional): Axes object to plot heatmap onto. If not specified, a new figure will be created. Defaults to None.
        guide_idx (Optional[Union[np.ndarray, List[bool]]], optional): Boolean array or list of indices to subset samples before calculating correlation. Defaults to None.
        correlation_prefix (str, optional): Prefix added to name of new columns added to AnnData object. Defaults to "".
        *args: Variable length argument list fed into sns.heatmap.
        **kwargs: Arbitrary keyword arguments fed into sns.heatmap.

    Returns:
        None
    """
    count_layer_label = count_layer or "X"
    corr_key = f"{correlation_prefix}corr_{count_layer_label}"
    if corr_key not in screen.obsm:
        _set_sample_correlation(
            screen, method, count_layer, guide_idx, correlation_prefix
        )
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        screen.obsm[corr_key],
        ax=ax,
        xticklabels=screen.obs.index,
        yticklabels=screen.obs.index,
        *args,
        **kwargs,
    )
    ax.set_box_aspect(1)
    ax.set_title(corr_key)


def sample_lfcs_correlation(
    screen: AnnData,
    cond1: str,
    cond2: str,
    rep_col: str = "replicate",
    cond_col: str = "condition",
    lognorm_counts_key: str = "lognorm_counts",
    method: str = "spearman",
    figsize: Tuple[int, int] = (6, 6),
    ax: Optional[matplotlib.axes.Axes] = None,
    guide_idx: Optional[Union[np.ndarray, List[bool]]] = None,
    *args,
    **kwargs,
) -> Tuple[pd.DataFrame, matplotlib.axes.Axes]:
    """Plot heatmap of pairwise correlation between guide log-fold-changes (LFCs) in response to two different conditions.

    Args:
        screen (AnnData): Annotated data matrix containing information about gene expression values.
        cond1 (str): Name of the first condition for comparison.
        cond2 (str): Name of the second condition for comparison.
        rep_col (str, optional): Column name for replicate identifier in screen.obs. Defaults to "rep".
        cond_col (str, optional): Column name for timepoint identifier in screen.obs. Defaults to "time".
        lognorm_counts_key (str, optional): Layer name for log-normalized counts in screen.layers. Defaults to "lognorm_counts".
        method (str, optional): Method used for calculating correlation. Defaults to "spearman".
        figsize ((int, int), optional): Width and height of figure in inches. Defaults to (10, 10).
        ax (matplotlib.axes.Axes, optional): Axes object to plot heatmap onto. If not specified, a new figure will be created. Defaults to None.
        guide_idx (Optional[Union[np.ndarray, List[bool]]], optional): Boolean array or list of indices to subset samples before calculating correlation. Defaults to None.
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        Tuple[pd.DataFrame, matplotlib.axes.Axes]: Pairwise correlation matrix between LFCs and axes object used for plotting.
    """
    if lognorm_counts_key not in screen.layers:
        raise ValueError(
            f"{lognorm_counts_key} not in screen.layers. Run pt.pp.log_norm(screen) to get log-normalized counts."
        )

    if guide_idx is not None:
        screen_subset = screen[:, guide_idx]
    else:
        screen_subset = screen
    lfcs = log_fold_change_reps(
        screen_subset,
        cond1=cond1,
        cond2=cond2,
        lognorm_counts_key=lognorm_counts_key,
        rep_col=rep_col,
        compare_col=cond_col,
    )

    if not ax:
        fig, ax = plt.subplots(figsize=figsize)
    if kwargs:
        if "annot" not in kwargs:
            kwargs["annot"] = True
    else:
        kwargs = {"annot": True}
    sns.heatmap(
        lfcs.corr(method=method).round(3),
        ax=ax,
        *args,
        **kwargs,
    )
    ax.set_box_aspect(1)
    return lfcs, ax
