
import pandas as pd
#import vintools as v
import re, numpy as np

def _transform_guide_spanning_dict_to_guide_df_addition(GuidesSpanning, guide_df):

    GuidePosition = {}
    for strand in ["+", "-"]:

        for protospacer, value in GuidesSpanning[strand].items():

            GuidePosition[protospacer] = [strand, value[0], value[1], value[2]]

    df_ = (
        pd.DataFrame(GuidePosition, index=["strand", "Chromosome", "Start", "End"])
        .T.reset_index()
        .rename({"index": "protospacer"}, axis=1)
    )
    
    guide_df = guide_df.merge(df_, on="protospacer", how="outer")
    guide_df["center"] = np.mean([guide_df["Start"], guide_df["End"]], axis=0)
    
    return guide_df

def _annotate_guide_position(GeneDict, guide_df):
        
    GuidesSpanning = v.ut.EmptyDict(["+", "-"])
    GuidesSpanning["+"] = {}
    GuidesSpanning["-"] = {}
    non_matching = []

    for target in GeneDict.keys():
        
        pos_target = GeneDict[target]["seq"]["+"]
        neg_target = GeneDict[target]["seq"]["-"]
        
        true_start = GeneDict[target]["Start"]*1e6
        true_end = GeneDict[target]["End"]*1e6

        target_guides = guide_df.loc[guide_df["target"] == target]
        
        for spacer in target_guides.protospacer.values:
            try:
                span = re.search(spacer, pos_target).span()
                guide_start, guide_end = span[0] + true_start, span[1] + true_start
                GuidesSpanning["+"][spacer] = (GeneDict[target]["Chromosome"], int(guide_start), int(guide_end))
            except:
                try:
                    span = re.search(spacer, neg_target).span()
                    guide_start, guide_end = true_end - span[0], true_end - span[1]
                    GuidesSpanning["-"][spacer] = (GeneDict[target]["Chromosome"], int(guide_start), int(guide_end))
                except:
                    non_matching.append([spacer, target])

    if len(non_matching) > 0:
        print("Number of non-mapping guides: {}".format(len(non_matching)))

    return _transform_guide_spanning_dict_to_guide_df_addition(GuidesSpanning, guide_df)