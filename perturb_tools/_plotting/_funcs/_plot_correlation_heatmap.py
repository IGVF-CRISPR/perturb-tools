
# _plot_correlation_heatmap.py
__module_name__ = "_plot_correlation_heatmap.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import os


# local imports #
# ------------- #
from ._get_default_matplotlib_figure_width_height import _get_default_matplotlib_figure_width_height
from ._set_matplotlib_rc_params import _set_matplotlib_rc_params


# setup figure params
_set_matplotlib_rc_params()


def _plot_correlation_heatmap(
    df, title="Sample Correlation", figsize=1, title_y=1.15, title_x=0, savename=False, cbar_aspect=20, cbar_shrink=0.6,
):

    default_h, default_w = _get_default_matplotlib_figure_width_height()
    
    fig = plt.figure(figsize=(default_w * figsize, default_h * figsize))
    gridspec = GridSpec(1, 1)
    ax = fig.add_subplot(gridspec[0, 0])
    im = ax.imshow(df.values.astype(float), cmap="plasma", vmin=0, vmax=1)
    plt.xticks(
        range(len(df.columns.values.astype(str))),
        df.columns.values.astype(str),
        rotation=60,
        ha="left",
    )
    ax.xaxis.tick_top()
    plt.yticks(
        range(len(df.columns.values.astype(str))),
        df.columns.values.astype(str),
        rotation=0,
    )
    plt.colorbar(im, shrink=cbar_shrink, aspect=cbar_aspect)
    plt.title(title, y=title_y, x=title_x, fontsize=16)

    if savename:
        figsavename = savename + ".png"
        plt.savefig(figsavename, bbox_inches="tight")
    plt.tight_layout()
    plt.show()
    if savename:
        return figsavename


