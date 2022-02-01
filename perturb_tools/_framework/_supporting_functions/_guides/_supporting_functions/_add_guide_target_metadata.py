
import pandas as pd

def _add_guide_target_metadata(guide_df, regex_annotations, DirectPairDict=False):

    """"""

    GuideMetadata = {}

    for guide_id in guide_df["barcode_id"].unique():
        if DirectPairDict:
            for key, value in DirectPairDict.items():
                if guide_id == key:
                    GuideMetadata[guide_id] = value
                else:
                    nullcount = 0
                    for annot in regex_annotations:
                        if annot in guide_id:
                            GuideMetadata[guide_id] = annot
                            break
                        else:
                            nullcount += 1
                            if nullcount == 3:
                                GuideMetadata[guide_id] = guide_id

    guide_metadata_df = (
        pd.DataFrame.from_dict(GuideMetadata, orient="index")
        .reset_index()
        .rename({"index": "barcode_id", 0: "target"}, axis=1)
    )

    guide_df = guide_df.merge(guide_metadata_df, on="barcode_id")

    return guide_df