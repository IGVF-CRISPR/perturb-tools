
#import vintools as v
import numpy as np
import matplotlib.pyplot as plt

def _cds_plot_construction():
    NotImplemented
    # plot = v.pl.ScatterPlot()
    # plot.construct_layout(figsize_width=4, figsize_height=0.5, nplots=1)
    # plot.style(grid=False)
    # ax = plot.AxesDict[0][0]

    # spines = v.pl.ax_spines(ax)
    # spines.delete(["left"])

    # return plot, ax


def _plot_cds(df, title, exon_spacing=100):

    """"""
    NotImplemented
    # c = v.pl.share_seq()["colors"]
    # c = np.repeat(c, 3).values

    # plot, ax = _cds_plot_construction()

    # floating_start = 0
    # exon_centers = []
    # for i in range(len(df)):
    #     stop = floating_start + df["exon_length"].iloc[i]
    #     exon_centers.append(np.mean([stop, floating_start]))
    #     ax.hlines(0.1, floating_start, stop, color=c[i], linewidth=10)
    #     floating_start = stop + exon_spacing

    # ax.set_ylabel(title, rotation=0, y=0.05)

    # ax.set_yticks([])
    # ax.set_ylim(0, 1)
    # xticks = ax.set_xticks(exon_centers)
    # xticks = ax.set_xticklabels(range(len(df)), fontsize=8)
    # plt.savefig(title + ".pdf")