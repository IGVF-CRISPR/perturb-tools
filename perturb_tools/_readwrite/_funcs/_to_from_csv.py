# _write_screen.py
__module_name__ = "_write_screen.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(
    [
        "vinyard@g.harvard.edu",
    ]
)


# package imports #
# --------------- #
import os
import pandas as pd
import anndata as ad

# local imports #
# ------------- #
from ..._utilities._funcs._flexible_mkdir import _flexible_multilevel_mkdir


def _save_dict_of_dfs(df_dict, out_path, filename, verbose=False):
    for key in df_dict.keys():
        file_outpath = os.path.join(out_path, "{}.{}.csv".format(filename, key))
        if verbose:
            print(file_outpath)
        df_dict[key].to_csv(file_outpath)


def _write_to_csv(screen, out_path="CRISPR_screen"):
    """This function will eventually be replaced with something more native / similar to AnnData."""

    _flexible_multilevel_mkdir(out_path)

    screen.X.tofile(
        os.path.join(out_path, "screen.X.csv"),
        sep=",",
    )
    screen.var.to_csv(os.path.join(out_path, "screen.var.csv"))
    screen.obs.to_csv(os.path.join(out_path, "screen.obs.csv"))
    _save_dict_of_dfs(screen.obsm, out_path, "screen.obsm")
    _save_dict_of_dfs(screen.obsp, out_path, "screen.obsp")
    _save_dict_of_dfs(screen.layers, out_path, "screen.layers")


def _read_from_csv(X_path=None, guide_path=None, sample_path=None, sep=","):
    if X_path is not None:
        X_df = pd.read_csv(X_path, delimiter=sep, header=0, index_col=0)
        X = X_df.values.T
    else:
        X = None
    guide_df = pd.read_csv(guide_path, sep=sep) if guide_path is not None else None
    sample_df = None if sample_path is None else pd.read_csv(sample_path, sep=sep)
    return ad.AnnData(X=X, var=guide_df, obs=sample_df)
