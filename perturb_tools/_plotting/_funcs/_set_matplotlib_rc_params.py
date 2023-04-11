import matplotlib
import os

def _set_matplotlib_rc_params():

    os.system("rm ~/.cache/matplotlib -rf")

    font = {"size": 12}
    matplotlib.rc(font)
    matplotlib.rcParams["font.sans-serif"] = "Arial"
    matplotlib.rcParams["font.family"] = "sans-serif"