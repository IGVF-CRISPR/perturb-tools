# _get_chromosome_sequence.py

from Bio import SeqIO


def _get_chromosome_sequence(ref_seq_path, query_chr, return_length=False):

    """
    Get a specific chromosome sequence from a reference genome. Also report the length of that sequence.
    
    Parameters:
    -----------
    ref_seq_path [ required ]
        Path to a reference genome fasta file.
    
    query_chr [ required ]
        Chromosome sequence to isolate
        type: str
        
    return_length [ optional ]
        default: False
        type: bool
        
    Returns:
    --------
    chromosome_reference_seq
        type: str
    
    len(chromosome_reference_seq) [ optional ]
    
    Notes:
    ------
    (1) if return_length is true, a list is returned to avoid setting two outputs.  
    """

    for record in SeqIO.parse(ref_seq_path, "fasta"):

        if record.description.split()[0] == query_chr:
            chromosome_reference_seq = str(record.seq)

        if return_length:
            return [chromosome_reference_seq, len(chromosome_reference_seq)]
        else:
            return chromosome_reference_seq

