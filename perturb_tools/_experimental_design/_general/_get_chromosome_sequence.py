
# _get_chromosome_sequence.py
__module_name__ = "_get_chromosome_sequence.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])

# package imports #
# --------------- #


def _get_chromosome_sequence(ref_seq_path, query_chr):

    """
    Get a specific chromosome sequence from a reference genome. Also report the length of that sequence.
    Parameters:
    -----------
    ref_seq_path
        Path to a reference genome fasta file.
    query_chr
        Chromosome sequence to isolate
        type: str

    Returns:
    --------
    chromosome_reference_seq
        type: str

    Notes:
    ------
    """

    for record in SeqIO.parse(ref_seq_path, "fasta"):

        if record.description.split()[0] == query_chr:
            chromosome_reference_seq = str(record.seq)

            return chromosome_reference_seq
