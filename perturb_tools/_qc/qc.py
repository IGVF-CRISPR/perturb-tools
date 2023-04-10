import numpy as np
import numpy.ma as ma
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
from perturb_tools import Screen


def set_count_stats(screen):
    screen.obs["median_X"] = np.nanmedian(screen.X, axis=1)
    screen.obs["median_X_bcmatch"] = np.nanmedian(screen.layers["X_bcmatch"], axis=1)
    if "edits" in screen.layers.keys():
        screen.obs["median_edit"] = np.nanmedian(screen.layers["edits"], axis=1)
        screen.get_guide_edit_rate()
        screen.obs["median_edit_rate"] = np.nanmedian(
            screen.layers["edit_rate"], axis=1
        )


def set_sample_correlation_guides(
    screen: Screen, guide_idx, prefix="", method="Pearson"
):
    corr_func_dict = {
        "Pearson": lambda X: ma.corrcoef(X.T),
        "Spearman": lambda X: scipy.stats.spearmanr(X, nan_policy="omit")[0],
    }
    corr_func = corr_func_dict[method]
    screen_subset = screen[guide_idx, :].copy()
    if prefix != "":
        prefix = f"{prefix}_"
    screen.obsm[f"{prefix}corr_X"] = corr_func(ma.masked_invalid(screen_subset.X))
    screen.obsm[f"{prefix}corr_X_bcmatch"] = corr_func(
        ma.masked_invalid(screen_subset.layers["X_bcmatch"])
    )
    screen.obsm[f"{prefix}corr_edits"] = corr_func(
        ma.masked_invalid(screen_subset.layers["edits"])
    )

    screen.obs[f"{prefix}mean_corr_X"] = screen.obsm[f"{prefix}corr_X"].mean(0)
    screen.obs[f"{prefix}mean_corr_X_bcmatch"] = screen.obsm[
        f"{prefix}corr_X_bcmatch"
    ].mean(0)
    screen.obs[f"{prefix}mean_corr_edits"] = screen.obsm[f"{prefix}corr_edits"].mean(0)

    screen.obs[f"{prefix}median_corr_X"] = np.median(
        screen.obsm[f"{prefix}corr_X"], axis=1
    )
    screen.obs[f"{prefix}median_corr_X_bcmatch"] = np.median(
        screen.obsm[f"{prefix}corr_X_bcmatch"], axis=1
    )
    screen.obs[f"{prefix}median_corr_edits"] = np.median(
        screen.obsm[f"{prefix}corr_edits"], axis=1
    )


def set_sample_correlation(screen: Screen, method="Pearson"):
    set_sample_correlation_guides(screen, screen.guides.index, prefix="", method=method)


def plot_correlation(screen, method="Pearson"):
    """
    arguments
    -- method: ["Pearson", "Spearman"]
    """
    set_sample_correlation(screen, method)
    n_corrs = len(screen.obsm.keys())
    fig, ax = plt.subplots(n_corrs, 1, figsize=(20, n_corrs * 20))
    for i, (k, cor) in enumerate(screen.obsm.items()):
        sns.heatmap(
            cor,
            ax=ax[i],
            xticklabels=screen.obs.index,
            yticklabels=screen.obs.index,
        )
        ax[i].set_title(k)


def G(v):
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


def set_X_gini(screen, plot=True):
    gini_vals = []
    if plot:
        plt.figure()
    bins = None
    for c in range(len(screen.obs.index)):
        bins, result, gini_val = G(screen.X[:, c])
        gini_vals.append(gini_val)
        if plot:
            plt.plot(
                bins,
                result,
                label="{} ({:.3f})".format(screen.obs.index[c], gini_val),
            )
    if plot:
        plt.plot(bins, bins, "--", label="perfect eq.")
        plt.xlabel("Cumulative fraction of guides \nfrom lowest to highest read counts")
        plt.ylabel("Cumulative fraction of read counts")
        plt.title("GINI index")
        plt.legend(bbox_to_anchor=(1.02, 1))
    screen.obs["gini_X"] = gini_vals


def set_X_bcmatch_gini(screen, plot=True):
    gini_vals = []
    if plot:
        plt.figure()
    bins = None
    for c in range(len(screen.obs.index)):
        bins, result, gini_val = G(screen.layers["X_bcmatch"][:, c])
        gini_vals.append(gini_val)
        if plot:
            plt.plot(
                bins,
                result,
                label="{} ({:.3f})".format(screen.obs.index[c], gini_val),
            )
    if plot:
        plt.plot(bins, bins, "--", label="perfect eq.")
        plt.xlabel("Cumulative fraction of guides \nfrom lowest to highest read counts")
        plt.ylabel("Cumulative fraction of read counts")
        plt.title("GINI index")
        plt.legend(bbox_to_anchor=(1.02, 1))
    screen.obs["gini_X_bcmatch"] = gini_vals
