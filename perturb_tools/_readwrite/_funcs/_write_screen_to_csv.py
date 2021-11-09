
# _write_screen.py
__module_name__ = "_write_screen.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import os


# local imports #
# ------------- #
from ..._utilities._funcs._flexible_mkdir import _flexible_multilevel_mkdir


def _save_dict_of_dfs(df_dict, out_path, filename, verbose=False):

    for key in df_dict.keys():
        file_outpath = os.path.join(out_path, "{}.{}.csv".format(filename, key))
        if verbose:
            print(file_outpath)
        df_dict[key].to_csv(file_outpath)


def _write_screen_to_csv(screen, out_path="CRISPR_screen"):

    """This function will eventually be replaced with something more native / similar to AnnData."""

    _flexible_multilevel_mkdir(out_path)

    screen.X.tofile(
        os.path.join(out_path, "screen.X.csv"),
        sep=",",
    )
    screen.guides.to_csv(os.path.join(out_path, "screen.guides.csv"))
    screen.condit.to_csv(os.path.join(out_path, "screen.condit.csv"))
    _save_dict_of_dfs(screen.condit_m, out_path, "screen.condit_m")
    _save_dict_of_dfs(screen.condit_p, out_path, "screen.condit_p")
    _save_dict_of_dfs(screen.layers, out_path, "screen.layers")
