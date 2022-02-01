class _SequenceManipulation:
    
    """
    Get the complement or reverse compliment of a DNA sequence.
    
    Parameters:
    -----------
    sequence
    
    Returns:
    --------
    complement_sequence
    
    reverse_sequence
    
    reverse_complement_sequence
    
    Notes:
    ------
    (1) No dependencies required. Pure python. 
    """
    
    def __init__(self, sequence):
        
        self.sequence = sequence
        self.ComplimentDict = {'C':'G', 'G':'C', 'T':'A', 'A':'T', 'N':'N'}
        self.complement_sequence = ""
        self.reverse_complement_sequence = ""
        
    def complement(self):
        
        """
        Get the complement of a DNA sequence.

        Parameters:
        -----------
        sequence (passed during instantiation).

        Returns:
        --------
        complement_sequence

        Notes:
        ------
        (1) No dependencies required. Pure python. 
        """
        for nucleotide in self.sequence:
            self.complement_sequence += self.ComplimentDict[nucleotide]
        
    
        return self.complement_sequence
        
    def reverse(self):
        
        """
        Get the reverse of a DNA sequence.

        Parameters:
        -----------
        sequence (passed during instantiation).

        Returns:
        --------
        reverse_sequence

        Notes:
        ------
        (1) No dependencies required. Pure python. 
        """
        
        self.reverse_sequence = self.sequence[::-1]
        
    def reverse_complement(self):
        
        """
        Get the reverse complement of a DNA sequence.

        Parameters:
        -----------
        sequence (passed during instantiation).

        Returns:
        --------
        reverse_complement_sequence

        Notes:
        ------
        (1) No dependencies required. Pure python. 
        """
        
        if not self.complement_sequence:
            for nucleotide in self.sequence:
                self.complement_sequence += self.ComplimentDict[nucleotide]
        
        self.reverse_complement_sequence = self.complement_sequence[::-1]
        return self.reverse_complement_sequence