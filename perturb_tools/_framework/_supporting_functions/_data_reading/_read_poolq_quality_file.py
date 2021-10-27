
import pandas as pd
import os

import vintools as v

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
    PoolQ_QualityDict = v.ut.EmptyDict(poolq_quality_keys)

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
    v.ut.mkdir_flex(analysis_dir)
    workbook_path = os.path.join(analysis_dir, workbook_path)
    v.ut.df_to_excel(
        list(quality_dfs.values()),
        workbook_path,
        sheetnames=poolq_quality_keys,
        index=False,
        silent=silent,
    )

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


def _read_file_in_dir(dir_path, filename):

    """Useful when you have multiple paths in a dir and you want to reuse the dir name. Uses the FileHandler class."""

    filepath = os.path.join(dir_path, filename)
    f = v.ut.FileHandler(filepath=filepath, verbose=False)
    return f.read(return_file=True)


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