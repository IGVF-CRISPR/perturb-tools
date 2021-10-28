import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# from ..._plotting._scatter._ScatterPlot_Module import ScatterPlot


def _calculate_LFC(df, condition, control):

    """
    
    df
        ResultsDict["lognorm_counts"]
        type: pandas DataFrame\
        
    condition        
    
    control
    
    """
    lfc = df[condition] - df[control]

    return lfc


def _get_conditions(lfc_df):
    return pd.Series(lfc_df.columns).str.split(" Rep", expand=True)[0]


def _calculate_multiLFC(df, conditions, controls):

    LFCDict = {}
    for i in range(len(conditions)):

        LFCDict[conditions[i] + " LFC"] = _calculate_LFC(
            df, condition=conditions[i], control=controls[i],
        )

    lfc_df = pd.DataFrame.from_dict(LFCDict)
    conditional_df = lfc_df.copy()
    conditional_df.columns = _get_conditions(conditional_df)

    return lfc_df, conditional_df


def _focus_guides(rowfile, guide_df, lfc_df_, condit_df_):

    row_df = pd.read_csv(
        rowfile, sep="\t", header=None, names=["sequence", "target_region"],
    )
    indexers = []

    for seq in guide_df.sequence:
        indexers.append(np.where(seq == row_df.sequence)[0][0])

    lfc_df = lfc_df_.iloc[np.array(indexers).astype(int)]
    condit_df = condit_df_.iloc[np.array(indexers).astype(int)]

    return lfc_df, condit_df


# def _plot_multiLFC(
#     df, rowfile, guide_df, conditions, controls, ms, elinewidth, narrow=True
# ):

#     lfc_df_, condit_df_ = _calculate_multiLFC(df, conditions, controls)

#     if narrow:
#         lfc_df, condit_df = _focus_guides(rowfile, guide_df, lfc_df_, condit_df_)
#     else:
#         lfc_df = lfc_df_
#         condit_df = condit_df_

#     LFC_Dict = {}

#     plot = ScatterPlot()
#     plot.construct_layout(
#         ncols=1,
#         figsize_width=20,
#         width_ratios=[1],
#         grid_hspace=0.25,
#         figsize_height=condit_df.columns.nunique(),
#         nplots=condit_df.columns.nunique(),
#     )
#     plot.style()

#     for n, cond in enumerate(condit_df.columns.unique()):
#         ax = plot.AxesDict[n][0]
#         LFC_Dict[cond] = {}
#         ax.set_title(cond, fontsize=50)
#         try:
#             LFC_Dict[cond]["LFC.mean"] = condit_df[cond].mean(axis=1)
#             LFC_Dict[cond]["LFC.stdev"] = condit_df[cond].std(axis=1)
#             ax.errorbar(
#                 guide_df.start,
#                 LFC_Dict[cond]["LFC.mean"],
#                 yerr=LFC_Dict[cond]["LFC.stdev"],
#                 fmt="o",
#                 ms=ms,
#                 elinewidth=elinewidth,
#             )

#         except:
#             ax.scatter(guide_df.start, condit_df[cond])

#     plt.show()
