import numpy as np
import numpy.ma as ma
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
from perturb_tools import Screen

def set_count_stats(screen):
    screen.condit['median_X'] = np.nanmedian(screen.X,axis=0)
    screen.condit['median_X_bcmatch'] = np.nanmedian(screen.layers['X_bcmatch'],axis=0)
    if 'edits' in screen.layers.keys():
        screen.condit['median_edit'] = np.nanmedian(screen.layers['edits'],axis=0)
        screen.get_guide_edit_rate()
        screen.condit['median_edit_rate'] = np.nanmedian(screen.layers['edit_rate'],axis=0)

def set_sample_correlation_guides(screen:Screen, guide_idx, prefix="", method="Pearson"):
    corr_func_dict = {
        "Pearson": lambda X: ma.corrcoef(X.T),
        "Spearman": lambda X: scipy.stats.spearmanr(X, nan_policy = 'omit')[0]
    }
    corr_func = corr_func_dict[method]
    screen_subset = screen[guide_idx,:].copy()
    if prefix != "": prefix = prefix + "_"
    screen.varm['{}corr_X'.format(prefix)] = corr_func(ma.masked_invalid(screen_subset.X))
    screen.varm['{}corr_X_bcmatch'.format(prefix)] = corr_func(ma.masked_invalid(screen_subset.layers["X_bcmatch"]))
    screen.varm['{}corr_edits'.format(prefix)] = corr_func(ma.masked_invalid(screen_subset.layers["edits"]))

    screen.condit['{}mean_corr_X'.format(prefix)] = screen.varm['{}corr_X'.format(prefix)].mean(0)
    screen.condit['{}mean_corr_X_bcmatch'.format(prefix)] = screen.varm['{}corr_X_bcmatch'.format(prefix)].mean(0)
    screen.condit['{}mean_corr_edits'.format(prefix)] = screen.varm['{}corr_edits'.format(prefix)].mean(0)

    screen.condit['{}median_corr_X'.format(prefix)] = np.median(screen.varm['{}corr_X'.format(prefix)], axis=0)
    screen.condit['{}median_corr_X_bcmatch'.format(prefix)] = np.median(screen.varm['{}corr_X_bcmatch'.format(prefix)], axis=0)
    screen.condit['{}median_corr_edits'.format(prefix)] = np.median(screen.varm['{}corr_edits'.format(prefix)], axis=0)

def set_sample_correlation(screen:Screen, method="Pearson"):
    set_sample_correlation_guides(screen, screen.guides.index, prefix="",method=method)

def plot_correlation(screen, method="Pearson"):
    """
    arguments
    -- method: ["Pearson", "Spearman"]
    """
    set_sample_correlation(screen, method)
    n_corrs = len(screen.varm.keys())
    fig, ax = plt.subplots(n_corrs, 1, figsize=(20, n_corrs * 20))
    for i, (k, cor) in enumerate(screen.varm.items()):
        sns.heatmap(cor, ax=ax[i], xticklabels=screen.condit.index, yticklabels=screen.condit.index)
        ax[i].set_title(k)

def G(v):
    bins = np.linspace(0., 1., 50)
    total = float(np.sum(v))
    yvals = []
    for b in bins:
        bin_vals = v[v <= np.quantile(v, b)]
        bin_fraction = (np.sum(bin_vals) / total)
        yvals.append(bin_fraction)
    # perfect equality area
    pe_area = np.trapz(bins, x=bins)
    # lorenz area
    lorenz_area = np.trapz(yvals, x=bins)
    gini_val = (pe_area - lorenz_area) / float(pe_area)
    return bins, yvals, gini_val

def set_X_gini(screen, plot = True):
    gini_vals = []
    if plot:
        plt.figure()
    for c in range(len(screen.condit.index)):
        bins, result, gini_val = G(screen.X[:,c])  
        gini_vals.append(gini_val)
        if plot:
            plt.plot(bins, result, label="{} ({:.3f})".format(screen.condit.index[c], gini_val))
    if plot:
        plt.plot(bins, bins, '--', label="perfect eq.")
        plt.xlabel("Cumulative fraction of guides \nfrom lowest to highest read counts")
        plt.ylabel("Cumulative fraction of read counts")
        plt.title("GINI index")
        plt.legend(bbox_to_anchor=(1.02, 1))
    screen.condit['gini_X'] = gini_vals
    

def set_X_bcmatch_gini(screen, plot = True):
    gini_vals = []
    if plot: plt.figure()

    for c in range(len(screen.condit.index)):
        bins, result, gini_val = G(screen.layers['X_bcmatch'][:,c])  
        gini_vals.append(gini_val)
        if plot: plt.plot(bins, result, label="{} ({:.3f})".format(screen.condit.index[c], gini_val))
    if plot:
        plt.plot(bins, bins, '--', label="perfect eq.")
        plt.xlabel("Cumulative fraction of guides \nfrom lowest to highest read counts")
        plt.ylabel("Cumulative fraction of read counts")
        plt.title("GINI index")
        plt.legend(bbox_to_anchor=(1.02, 1))
    screen.condit["gini_X_bcmatch"] = gini_vals
    