# _read_screen_from_PoolQ.py
__module_name__ = "_read_screen_from_PoolQ.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

import pandas as pd
import numpy as np
import vintools as v

from ._read_poolq_quality_file import _read_poolq_quality_file

def _assemble_pandas_dict(filepaths, keys=False, sep="\t", **kwargs):

    """"""

    PandasDict = {}
    for i, path in enumerate(filepaths):
        if keys:
            df_key = keys[i]
        else:
            df_key = keys[i]
        PandasDict[df_key] = pd.read_csv(path, sep, **kwargs)

    return PandasDict


def _find_key_path(outs_dict, key):
    
    key_loc = np.where(np.array(list(outs_dict.keys())) == key)[0][0] - 1

    return outs_dict["paths"][key_loc]


def _read_poolq_counts_dfs(
    poolq_outs, counts_keys=["counts", "lognormalized-counts"], **kwargs
):

    """"""

    filepaths = []
    for key in counts_keys:
        filepaths.append(_find_key_path(poolq_outs, key))
    
    pd_dict = _assemble_pandas_dict(filepaths, keys=counts_keys, sep="\t", **kwargs)

    return pd_dict

def _load_parse_PoolQ_counts_df(
    PoolQ_OutsDict, ScreenDict=False, counts_keys=["counts", "lognormalized-counts"]
):

    """
    Load counts df (raw and lognorm) from PoolQ Dict. 
    
    Parameters:
    -----------
    PoolQ_OutsDict
    
    ScreenDict
    
    counts_keys
    
    Returns:
    --------
    ScreenDict
    
    Notes:
    ------
    """

    if not ScreenDict: # initialize the dict if it does not exist
        ScreenDict = v.ut.EmptyDict(["layers", "condit_p", "condit_m", "uns"])
        
    df_dict = _read_poolq_counts_dfs(PoolQ_OutsDict, counts_keys)

    ScreenDict["X"] = df_dict["counts"].drop(["Row Barcode", "Row Barcode IDs"], axis=1).values
    ScreenDict["guides"] = df_dict["counts"][['Row Barcode', 'Row Barcode IDs']]
    ScreenDict["guides"].columns = ['barcode', 'barcode_id']
    ScreenDict["layers"]["lognorm_counts"] = df_dict["lognormalized-counts"].drop(["Row Barcode", "Row Barcode IDs"], axis=1)
    ScreenDict["condit"] = pd.DataFrame(ScreenDict['layers']["lognorm_counts"].columns, columns=["Condition"])

    return ScreenDict

def _compose_condit_mp(
    ScreenDict, PoolQ_OutsDict, counts_keys=["barcode-counts", "unexpected-sequences"]
):

    """
    Load counts df (raw and lognorm) from PoolQ Dict.

    Parameters:
    -----------
    PoolQ_OutsDict

    ScreenDict

    counts_keys

    Returns:
    --------
    ScreenDict

    Notes:
    ------
    """
    
    ScreenDict["condit_m"] = _read_poolq_counts_dfs(
        PoolQ_OutsDict, counts_keys, index_col=[0, 1, -1]
    )
    ScreenDict["condit_p"] = _read_poolq_counts_dfs(PoolQ_OutsDict, ["correlation"], index_col=0)

    return ScreenDict


def _compose_uns_from_quality_file(ScreenDict, PoolQ_OutsDict, poolq_outs_path, out_path):

    quality_dfs = _read_poolq_quality_file(
        poolq_outs_path, analysis_dir=out_path, return_df=True
    )
    
    ScreenDict["uns"]["run_info"] = PoolQ_OutsDict["runinfo"]
    ScreenDict["uns"]["poolq3"] = PoolQ_OutsDict["poolq3"]

    for key, _df in quality_dfs.items():
        ScreenDict["uns"][str(key)] = _df

    return ScreenDict

def _read_screen_from_PoolQ(poolq_outs_path, analaysis_out_path="wokrbook.xlsx"):

    PoolQ_OutsDict = v.ut.glob_dict(poolq_outs_path)
    
    ScreenDict = _load_parse_PoolQ_counts_df(PoolQ_OutsDict, ScreenDict=False, counts_keys=["counts", "lognormalized-counts"])
    ScreenDict = _compose_condit_mp(ScreenDict, PoolQ_OutsDict, counts_keys=["barcode-counts", "unexpected-sequences"])
    ScreenDict = _compose_uns_from_quality_file(ScreenDict, PoolQ_OutsDict, poolq_outs_path, out_path=analaysis_out_path)
    
    return ScreenDict
