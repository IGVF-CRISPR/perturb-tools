
from ....._utilities._funcs._SequenceManipulation import _SequenceManipulation
from ....._utilities._funcs._fetch_chromosome_sequence import _fetch_chromosome_sequence

def _load_relevant_chromosomes(gene_df, ref_seq_path, silent=False):

    """
    Load unique chromosomes found in the gene dict into memory


    Notes:
    ------
    (1) Doing this up front and once prevents the redundant fetching of chromosome sequences, which can be slow.
    (2) The produced dict can be deleted (and the corresponding memory freed) upon substringing the chromosome.
    """

    ChromDict = {}

    for chrom in gene_df.Chromosome.unique():
        if not silent:
            print("loading {}...".format(chrom))
        ChromDict[chrom] = _fetch_chromosome_sequence(ref_seq_path, chrom)
    return ChromDict

def _add_context_sequence_to_gene_dict(gene_dict, gene_df, ref_seq_path):

    ChromDict = _load_relevant_chromosomes(gene_df, ref_seq_path)

    for gene, values in gene_dict.items():
        gene_dict[gene]["seq"] = {}
        gene_dict[gene]["chrom.length"] = len(ChromDict[values["Chromosome"]])
        gene_start = int(values["Start"] * 1e6)
        gene_end = int(values["End"] * 1e6)
        sequence = _SequenceManipulation(
            ChromDict[values["Chromosome"]][gene_start:gene_end]
        )

        gene_dict[gene]["seq"]["+"] = sequence.sequence
        gene_dict[gene]["seq"]["-"] = sequence.reverse_complement()

    print("\nDeleteting whole chromosome sequences from memory...")
    del ChromDict

    return gene_dict