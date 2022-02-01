
from ..._arithmetic._funcs._SequenceManipulation import _SequenceManipulation
from ..._utilities._funcs._scan_sequence_for_motif import _scan_sequence_for_motif

import numpy as np

def _scan_for_BsmbI(df, filter_BsmbI=True):

    """"""

    forward_motif_search = "TCTC"
    reverse_motif_search = "GAGAC"

    BsmbI_flag = []
    count = 0
    for a, guide in enumerate(df.guide_context):
        
        
        guide_seq = _SequenceManipulation(guide)
        rev_guide = guide_seq.reverse_complement()

        a = _scan_sequence_for_motif(guide, forward_motif_search, negative_report=True)
        b = _scan_sequence_for_motif(guide, reverse_motif_search, negative_report=True)
        c = _scan_sequence_for_motif(
            rev_guide, forward_motif_search, negative_report=True
        )
        d = _scan_sequence_for_motif(
            rev_guide, reverse_motif_search, negative_report=True
        )

        x = np.array([a, b, c, d])

        if np.any(x) == True:
            count += 1
            BsmbI_flag.append(True)
        else:
            BsmbI_flag.append(False)
    print(count, "sgRNAs were found to contain BsmbI sites and were flagged.")

    df["BsmbI_flag"] = BsmbI_flag
    
    if filter_BsmbI:
        print("\nFiltering BsmbI-containing sgRNAs...\n")
        BsmbI_guide_df = df.loc[df.BsmbI_flag == True]
        df = df.loc[df.BsmbI_flag == False]
        return df, BsmbI_guide_df

    return df, None