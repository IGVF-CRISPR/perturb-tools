
def _scan_sequence_for_motif(seq, motif, negative_report=False):

    """Return location of motif if found"""
    
    if motif in seq:
        return seq.find(motif)
    if negative_report == True:
        return False