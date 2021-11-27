# _read_screen_from_PoolQ.py

__module_name__ = "_read_screen_from_PoolQ.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# def _assemble_pandas_dict(filepaths, keys=False, sep="\t", **kwargs):

#     """"""

#     PandasDict = {}
#     for i, path in enumerate(filepaths):
#         if keys:
#             df_key = keys[i]
#         else:
#             df_key = keys[i]
#         PandasDict[df_key] = pd.read_csv(path, sep, **kwargs)

#     return PandasDict


# def _find_key_path(outs_dict, key):
    
#     key_loc = np.where(np.array(list(outs_dict.keys())) == key)[0][0] - 1

#     return outs_dict["paths"][key_loc]


# def _read_poolq_counts_dfs(
#     poolq_outs, counts_keys=["counts", "lognormalized-counts"], **kwargs
# ):

#     """"""

#     filepaths = []
#     for key in counts_keys:
#         filepaths.append(_find_key_path(poolq_outs, key))
    
#     pd_dict = _assemble_pandas_dict(filepaths, keys=counts_keys, sep="\t", **kwargs)

#     return pd_dict

# def _load_parse_PoolQ_counts_df(
#     PoolQ_OutsDict, ScreenDict=False, counts_keys=["counts", "lognormalized-counts"]
# ):

#     """
#     Load counts df (raw and lognorm) from PoolQ Dict. 
    
#     Parameters:
#     -----------
#     PoolQ_OutsDict
    
#     ScreenDict
    
#     counts_keys
    
#     Returns:
#     --------
#     ScreenDict
    
#     Notes:
#     ------
#     """

#     if not ScreenDict: # initialize the dict if it does not exist
#         ScreenDict = v.ut.EmptyDict(["layers", "condit_p", "condit_m", "uns"])
        
#     df_dict = _read_poolq_counts_dfs(PoolQ_OutsDict, counts_keys)

#     ScreenDict["X"] = df_dict["counts"].drop(["Row Barcode", "Row Barcode IDs"], axis=1).values
#     ScreenDict["guides"] = df_dict["counts"][['Row Barcode', 'Row Barcode IDs']]
#     ScreenDict["guides"].columns = ['barcode', 'barcode_id']
#     ScreenDict["layers"]["lognorm_counts"] = df_dict["lognormalized-counts"].drop(["Row Barcode", "Row Barcode IDs"], axis=1)
#     ScreenDict["condit"] = pd.DataFrame(ScreenDict['layers']["lognorm_counts"].columns, columns=["Condition"])

#     return ScreenDict

# def _compose_condit_mp(
#     ScreenDict, PoolQ_OutsDict, counts_keys=["barcode-counts", "unexpected-sequences"]
# ):

#     """
#     Load counts df (raw and lognorm) from PoolQ Dict.

#     Parameters:
#     -----------
#     PoolQ_OutsDict

#     ScreenDict

#     counts_keys

#     Returns:
#     --------
#     ScreenDict

#     Notes:
#     ------
#     """
    
#     ScreenDict["condit_m"] = _read_poolq_counts_dfs(
#         PoolQ_OutsDict, counts_keys, index_col=[0, 1, -1]
#     )
#     ScreenDict["condit_p"] = _read_poolq_counts_dfs(PoolQ_OutsDict, ["correlation"], index_col=0)

#     return ScreenDict


# def _compose_uns_from_quality_file(ScreenDict, PoolQ_OutsDict, poolq_outs_path, out_path):

#     quality_dfs = _read_poolq_quality_file(
#         poolq_outs_path, analysis_dir=out_path, return_df=True
#     )
    
#     ScreenDict["uns"]["run_info"] = PoolQ_OutsDict["runinfo"]
#     ScreenDict["uns"]["poolq3"] = PoolQ_OutsDict["poolq3"]

#     for key, _df in quality_dfs.items():
#         ScreenDict["uns"][str(key)] = _df

#     return ScreenDict

# def _read_screen_from_PoolQ(poolq_outs_path, analaysis_out_path="wokrbook.xlsx"):

#     PoolQ_OutsDict = v.ut.glob_dict(poolq_outs_path)
    
#     ScreenDict = _load_parse_PoolQ_counts_df(PoolQ_OutsDict, ScreenDict=False, counts_keys=["counts", "lognormalized-counts"])
#     ScreenDict = _compose_condit_mp(ScreenDict, PoolQ_OutsDict, counts_keys=["barcode-counts", "unexpected-sequences"])
#     ScreenDict = _compose_uns_from_quality_file(ScreenDict, PoolQ_OutsDict, poolq_outs_path, out_path=analaysis_out_path)
    
#     return ScreenDict

##################

import pandas as pd
import numpy as np
import os

from ..._utilities._funcs._glob_dict import _glob_dict
from ..._utilities._funcs._FileHandler import _FileHandler

def _read_file_in_dir(dir_path, filename):

    """Useful when you have multiple paths in a dir and you want to reuse the dir name. Uses the FileHandler class."""

    filepath = os.path.join(dir_path, filename)
    f = _FileHandler(filepath=filepath, verbose=False)
    return f.read(return_file=True)

def _create_EmptyDict(keys):

    """Create an empty python Dictionary given a set of keys."""

    Dict = {}
    for key in keys:
        Dict[key] = {}

    return Dict

def _read_poolq_quality_file(
    poolq_outs_path,
    analysis_dir,
    workbook_path="quality.xlsx",
    return_df=False,
    quality_filename="quality.txt",
    poolq_quality_keys=[
        "metadata",
        "bc_readcounts",
        "common_barcodes",
    ],
    verbose=False,
    silent=True,
):

    """"""

    QualityDataFrames = {}

    quality_file = _read_file_in_dir(poolq_outs_path, quality_filename)
    PoolQ_QualityDict = _create_EmptyDict(poolq_quality_keys)

    for n, line in enumerate(quality_file):
        line = _parse_line(line)
        if not line is None:
            if ":" in line[0]:
                PoolQ_QualityDict["metadata"][n] = line[0].split(":")
            elif len(line) > 2:
                _read_count_for_df(
                    PoolQ_QualityDict["bc_readcounts"], n, line
                )
            elif len(line) == 2:
                _read_count_for_df(
                    PoolQ_QualityDict["common_barcodes"], n, line
                )
            else:
                if verbose:
                    print("WARNING: LINE NOT ASSIGNED: {}".format(line))

    PoolQ_QualityDict["metadata"]["cols"] = ["metric", "statistic"]
    PoolQ_QualityDict["common_barcodes"]["cols"] = ["barcode", "count"]

    quality_dfs = {}
    for key in poolq_quality_keys:
        quality_dfs[key] = _make_df_with_colnames(PoolQ_QualityDict[key])
#     v.ut.mkdir_flex(analysis_dir)
    workbook_path = os.path.join(analysis_dir, workbook_path)
#     v.ut.df_to_excel(
#         list(quality_dfs.values()),
#         workbook_path,
#         sheetnames=poolq_quality_keys,
#         index=False,
#         silent=silent,
#     )

    if return_df:
        return quality_dfs
    
    
def _assemble_PoolQ_Quality_Dict():

    return _create_EmptyDict(
        keys=["metadata", "bc_readcounts", "common_barcodes"]
    )


def _make_df_with_colnames(dictionary, colnames_key="cols"):

    return (
        pd.DataFrame.from_dict(
            dictionary, orient="index", columns=dictionary[colnames_key]
        )
        .drop(colnames_key)
        .reset_index(drop=True)
    )


def _parse_line(line):
    """Parse poolq Quality file lines."""
    if not line == "\n":
        parsed_line = line.strip("\n").split("\t")
        if not parsed_line[0].startswith("Read counts"):
            return parsed_line


def _read_count_for_df(CountDict, n, line):

    """"""

    if len(CountDict) == 0:
        CountDict["cols"] = line
    else:
        CountDict[n] = line

    return CountDict

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
        ScreenDict = _create_EmptyDict(["layers", "condit_p", "condit_m", "uns"])
        
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

    PoolQ_OutsDict = _glob_dict(poolq_outs_path)
    
    ScreenDict = _load_parse_PoolQ_counts_df(PoolQ_OutsDict, ScreenDict=False, counts_keys=["counts", "lognormalized-counts"])
    ScreenDict = _compose_condit_mp(ScreenDict, PoolQ_OutsDict, counts_keys=["barcode-counts", "unexpected-sequences"])
    ScreenDict = _compose_uns_from_quality_file(ScreenDict, PoolQ_OutsDict, poolq_outs_path, out_path=analaysis_out_path)
    
    return ScreenDict


