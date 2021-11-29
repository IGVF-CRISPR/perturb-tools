
import numpy as np
import pandas as pd
from ..._arithmetic._funcs._SequenceManipulation import _SequenceManipulation

def _add_protospacer_sequence(df, chromosome_seq, strand):

    """Add the sequence of the protospacer and the PAM to the df."""

    df["PAM_loci"] = df.PAM_loci.astype(int)

    pam_seqs = np.array([])
    guide_seqs = np.array([])

    if strand == "+":

        for i in df['PAM_loci']:
            pam_seqs = np.append(pam_seqs, chromosome_seq[i : i + 3])
        for i in df.index:
            guide_seqs = np.append(
                guide_seqs, chromosome_seq[df.guide_begin[i] : df.guide_end[i]]
            )

    elif strand == "-":

        for i in df['PAM_loci']:            
            pam_seq = _SequenceManipulation(chromosome_seq[i - 3 : i])
            pam_seqs = np.append(pam_seqs, pam_seq.reverse_complement())
                
        for i in df.index:
            guide_seq = _SequenceManipulation(chromosome_seq[df.guide_end[i] : df.guide_begin[i]])
            guide_seqs = np.append(guide_seqs, guide_seq.reverse_complement())
    else:
        print("Strand unassigned. Interrupting.")
        return

    df["PAM"] = pam_seqs
    df["protospacer"] = guide_seqs

    return df

def _annotate_protospacer(df, chromosome_seq):

    """
    Annotate protospacers and PAM sequence. 
    
    Annotate df with information about guide protospacer loci and sequence.
    Splits df into forward and reverse strand, adds loci of protospacer then recombines as before.
    """
    forw_df = df.loc[df.strand == "+"]
    revr_df = df.loc[df.strand == "-"]

    forw_df["guide_begin"] = (forw_df.PAM_loci - 20).astype(int)
    forw_df["guide_end"] = (forw_df.PAM_loci).astype(int)
    forw_df = _add_protospacer_sequence(forw_df, chromosome_seq, strand="+")

    revr_df["guide_begin"] = (revr_df.PAM_loci + 20).astype(int)
    revr_df["guide_end"] = (revr_df.PAM_loci).astype(int)
    revr_df = _add_protospacer_sequence(revr_df, chromosome_seq, strand="-")

    df = pd.concat([forw_df, revr_df]).sort_values("PAM_loci").reset_index(drop=True)

    return df