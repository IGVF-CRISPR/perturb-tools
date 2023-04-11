import numpy as np
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/39512260/calculating-gini-coefficient-in-python-numpy
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

def plot_gini(adata, ax):
    for c in adata.condit.index:
        bins, result, gini_val = G(adata[:,c].X)  
        ax.plot(bins, result, label="{} ({:.3f})".format(c, gini_val))
    ax.plot((0, 1), (0, 1), '--', label="perfect eq.")
    ax.set_xlabel("Cumulative fraction of guides \nfrom lowest to highest read counts")
    ax.set_ylabel("Cumulative fraction of read counts")
    ax.set_title("GINI index")
    ax.set_legend(bbox_to_anchor=(1.02, 1))


def plot_count(adata, ax):
    ax.figure()
    bins = np.linspace(0, adata.X.max(), 100)
    for c in adata.condit.index:
        ax.hist(adata[:,c].X, bins, alpha = 0.5, 
                label="{} (median={})".format(c, np.nanmedian(adata[:,c].X)))
    ax.set_xlabel("Cumulative fraction of guides \nfrom lowest to highest read counts")
    ax.set_ylabel("Cumulative fraction of read counts")
    ax.set_title("GINI index")
    ax.set_legend(bbox_to_anchor=(1.02, 1))