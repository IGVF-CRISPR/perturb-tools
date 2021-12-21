# _fetch_chromosome_sequence.py

from Bio import SeqIO

def _fetch_chromosome_sequence(ref_seq_path, query_chr):
    
    """
    Isolate a specific chromosomal sequence (forward annotation) from a reference genome. 
    
    ref_seq_path
        Path to reference genome sequence.
        type: str
        
    query_chr
        example: 'chr2'
        type: str
        
    silent
        default: False
        default: boolean
    """
    
    for record in SeqIO.parse(ref_seq_path, "fasta"):
        if record.description.split()[0] == query_chr:
            return str(record.seq)
