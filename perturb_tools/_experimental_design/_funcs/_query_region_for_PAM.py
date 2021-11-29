
from ..._utilities._funcs._SequenceManipulation import _SequenceManipulation
from Bio import SeqIO
import pandas as pd
import numpy as np
import regex


def _get_seq_one_chr(ref_seq_path, query_chr):

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
    
    len(chromosome_reference_seq)
    """

    for record in SeqIO.parse(ref_seq_path, "fasta"):

        if record.description.split()[0] == query_chr:
            chromosome_reference_seq = str(record.seq)

            return chromosome_reference_seq, len(chromosome_reference_seq)

def _adjust_reverse_strand_coords_to_forward(
    reverse_seq_pam_match_coords, chromosome_length
):

    """
    Given a pam and a sequece, find sgRNAs within that sequence.
    Parameters:
    -----------
    reverse_seq_pam_match_coords
        Numpy array contianing loci of the first base in the PAM sequence. Numbers start right to left (5' to 3')
    Returns:
    --------
    reverse_seq_pam_match_coords_adjusted
        Numpy array contianing loci of the first base in the PAM sequence. Numbers start left to right (3' to 5')
    """

    return (chromosome_length - reverse_seq_pam_match_coords)


def _find_pam_loci_whole_chromosome(ref_seq_path, query_chr, PAM="NGG"):

    """
    Given a pam and a sequece, find sgRNAs within that sequence.
    Parameters:
    -----------
    seq: DNA sequence
    PAM: PAM sequence. Default=='NGG'
    Returns:
    --------
    forward_PAM_loci
        Numpy array contianing loci of the first base in the PAM sequence. Numbers start left to right (5' to 3')
    reverse_PAM_loci
        Numpy array contianing loci of the first base in the PAM sequence. Numbers start right to left (5' to 3')
    """

    reference_chromosome_sequence, chromosome_length = _get_seq_one_chr(
        ref_seq_path, query_chr=query_chr
    )
    print(
        "Searching",
        str(chromosome_length),
        "bp along",
        query_chr,
        "for PAM sequences...",
    )

    pam = PAM.replace("N", ".")

    forward_seq_pam_match_coords = [
        m.span()[0]
        for m in regex.finditer(
            r"%s" % (pam), reference_chromosome_sequence, overlapped=True
        )
    ]
    
    ref_seq = _SequenceManipulation(reference_chromosome_sequence)
    reverse_complement_seq = ref_seq.reverse_complement()
    
    reverse_complement_seq_pam_match_coords = [
        m.span()[0]
        for m in regex.finditer(r"%s" % (pam), reverse_complement_seq, overlapped=True)
    ]

    forward_seq_pam_match_coords = np.array(forward_seq_pam_match_coords)
    reverse_complement_seq_pam_match_coords = np.array(
        reverse_complement_seq_pam_match_coords
    )

    print(
        str(len(forward_seq_pam_match_coords)),
        "\nPAM matches identified along the forward strand of query chromsome.\n",
    )
    print(
        str(len(reverse_complement_seq_pam_match_coords)),
        "\nPAM matches identified along the reverse strand of query chromsome.\n",
    )
    
    reverse_complement_seq_pam_match_coords_adjusted = (
        _adjust_reverse_strand_coords_to_forward(
            reverse_complement_seq_pam_match_coords, chromosome_length
        )
    )

    return (
        forward_seq_pam_match_coords,
        reverse_complement_seq_pam_match_coords_adjusted,
        reference_chromosome_sequence
    )

def _get_PAMS_in_regions(regions_df, ref_seq_path, chromosome_key="Chromosome", PAM="NGG"):
    
    """
    Looks for PAM loci one chromosome at a time. Any chromosome specified in the peaks bed file is included.
    
    Parameters:
    -----------
    atac_peaks_path
        Path to a bed file of genomic loci. 
    
    ref_seq_path
        Path to a reference genome fasta file.
    PAM
        PAM motif to search (default: NGG)
        type: str
    Returns:
    --------
    """

    forward_pams = []
    reverse_pams = []
    chrom_seqs = []

    for chrom in regions_df[chromosome_key].unique():

        forw_chr, revr_chr, chrom_seq = _find_pam_loci_whole_chromosome(
            ref_seq_path, query_chr=chrom, PAM=PAM
        )
        forward_pams.append(forw_chr)
        reverse_pams.append(revr_chr)
        chrom_seqs.append(chrom_seq)

    return forward_pams, reverse_pams, chrom_seqs

def _query_region_for_PAM(regions_df, ref_seq_path, chromosome_key="Chromosome", PAM="NGG"):

    """
    Looks for PAM loci one chromosome at a time. Any chromosome specified in the peaks bed file is included.

    Parameters:
    -----------
    regions_df
        Path to a bed file of genomic loci.

    ref_seq_path
        Path to a reference genome fasta file.
    PAM
        PAM motif to search (default: NGG)
        type: str
    Returns:
    --------


    """

    forward_pams = []
    reverse_pams = []
    chrom_seqs = []

    for chrom in regions_df[chromosome_key].unique():

        forw_chr, revr_chr, chrom_seq = _find_pam_loci_whole_chromosome(
            ref_seq_path, query_chr=chrom, PAM=PAM
        )
        forward_pams.append(forw_chr)
        reverse_pams.append(revr_chr)
        chrom_seqs.append(chrom_seq)
        
    return forward_pams, reverse_pams, chrom_seqs