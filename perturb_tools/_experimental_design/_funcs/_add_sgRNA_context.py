### get context info for sgRNA

import pandas as pd

def _get_reverse_sgRNA_context(seq, df, flank_N):

    """"""

    # isolate all reverse strand sgRNAs
    reverse_guides = df.loc[df.strand == "-"]

    context_bound_5prime = reverse_guides.PAM_loci - [3 + flank_N]
    context_bound_3prime = reverse_guides.PAM_loci + [20 + flank_N]

    guide_context_seqs = []

    if len(context_bound_5prime) == len(context_bound_3prime):
        for i in range(len(context_bound_5prime)):
            guide_context = seq[
                context_bound_5prime.tolist()[i] : context_bound_3prime.tolist()[i]
            ]
            guide_context_seqs.append(guide_context)

    reverse_guides["guide_context"] = guide_context_seqs

    return reverse_guides


def _get_forward_sgRNA_context(seq, df, flank_N):

    """"""

    # isolate all reverse strand sgRNAs
    forward_guides = df.loc[df.strand == "+"]

    context_bound_5prime = forward_guides.PAM_loci - [20 + flank_N]
    context_bound_3prime = forward_guides.PAM_loci + [3 + flank_N]

    guide_context_seqs = []

    if len(context_bound_5prime) == len(context_bound_3prime):
        for i in range(len(context_bound_5prime)):
            guide_context = seq[
                context_bound_5prime.tolist()[i] : context_bound_3prime.tolist()[i]
            ]
            guide_context_seqs.append(guide_context)

    forward_guides["guide_context"] = guide_context_seqs

    return forward_guides


def _add_sgRNA_context(df, chrom_seq, flank_N=10):

    """
    Get the flanking N nts (default N=10) around a protospacer + PAM
    DF and chrom_seqs must be sorted by chromosome order
    """

    dfs = []

    for i, chrom in enumerate(df["Chromosome"].unique()):

        chrom_df = df.loc[df["Chromosome"] == chrom]

        forw_context_added_df = _get_forward_sgRNA_context(
            chrom_seq[i], chrom_df, flank_N
        )
        rev_context_added_df = _get_reverse_sgRNA_context(
            chrom_seq[i], chrom_df, flank_N
        )

        dfs.append(forw_context_added_df)
        dfs.append(rev_context_added_df)

    df = pd.concat(dfs)

    return df