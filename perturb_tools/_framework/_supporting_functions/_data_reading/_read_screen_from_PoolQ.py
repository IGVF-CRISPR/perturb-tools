
# _read_screen_from_PoolQ.py
__module_name__ = "_read_screen_from_PoolQ.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])



# def _read_poolq(poolq_outs_path):

#     """

#     Returns:
#     --------
#     screen

#     Notes:
#     ------
#     (1) screen is an "AnnData like" option.

#     """

#     files = glob.glob(poolq_outs_path + "*")

#     AllFiles = {}
#     AllFiles["full_paths"] = []

#     for filepath in files:
#         AllFiles["full_paths"].append(filepath)
#         f_ = FileHandler(filepath)
#         AllFiles[f] = f_.read(return_file=True)

#     return AllFiles


def _read_screen_from_PoolQ(poolq_outs_path):

    PoolQ_outs = _glob_dict(poolq_outs_path)

    counts_df = pd.read_csv(
        PoolQ_outs["paths"][4], sep="\t", index_col=[0, 1]
    )  # has multiindex

    screen = {}
    screen["X"] = counts_df.values
    screen["guides"] = counts_df.index.to_frame().reset_index(drop=True)
    screen["condit"] = pd.DataFrame(counts_df.columns[2:], columns=["conditions"])
    screen["layers"] = {}
    screen["layers"]["X_lognorm"] = pd.read_csv(
        PoolQ_outs["paths"][3], sep="\t", index_col=[0, 1]
    )

    screen["condit_p"] = {}
    screen["condit_m"] = {}
    screen["uns"] = {}

    screen["condit_m"]["barcode_counts"] = pd.read_csv(
        PoolQ_outs["paths"][1], sep="\t", index_col=[0, 1, -1]
    )  # top 100 only?
    screen["condit_m"]["unexpected_sequences"] = pd.read_csv(
        PoolQ_outs["paths"][2], sep="\t", index_col=[0, 1, -1]
    )

    quality_dfs = _read_poolq_quality_file(
        poolq_outs_path, analysis_dir="./newouts", return_df=True
    )
    screen["condit_p"]["correlation"] = pd.read_csv(
        PoolQ_outs["paths"][5], sep="\t", index_col=0
    )

    screen["uns"]["run_info"] = PoolQ_outs["runinfo"]
    screen["uns"]["poolq3"] = PoolQ_outs["poolq3"]
    for key, _df in quality_dfs.items():
        screen["uns"][str(key)] = _df

    return screen
