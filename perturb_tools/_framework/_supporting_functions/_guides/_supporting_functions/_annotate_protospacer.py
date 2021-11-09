import numpy as np

def _annotate_protospacer(guide_df, protospacer_length=20):

    if np.all(guide_df["barcode"].str.len().unique() == protospacer_length):
        guide_df["protospacer"] = guide_df["barcode"]

    return guide_df
