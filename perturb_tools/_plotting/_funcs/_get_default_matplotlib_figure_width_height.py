
# _get_default_matplotlib_figure_width_height.py
__module_name__ = "_get_default_matplotlib_figure_width_height.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import matplotlib, numpy


def _get_default_matplotlib_figure_width_height():

    """
    Return default height and width of matplotlib figures.
    
    Parameters:
    -----------
    None
    
    Returns:
    --------
    w, h
        width and height of the default matplotlib figure size
    
    Notes:
    ------
    
    """

    default_wh = matplotlib.rcParams["figure.figsize"]  # w x h
    w, h = default_wh[0], default_wh[1]

    return numpy.array(w, h)