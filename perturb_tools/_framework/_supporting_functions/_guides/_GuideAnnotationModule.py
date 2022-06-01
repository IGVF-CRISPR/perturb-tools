
import pandas as pd

from ._supporting_functions._annotate_protospacer import _annotate_protospacer
from ._supporting_functions._add_guide_target_metadata import _add_guide_target_metadata
from ._supporting_functions._add_context_sequence_to_gene_dict import _add_context_sequence_to_gene_dict

class _GuideAnnotation:
    def __init__(self, guide_df, genes, chromosome, start, stop):

        """"""

        self.guide_df = _annotate_protospacer(guide_df)
        self.chr = chromosome
        self.start = start
        self.stop = stop
        #self.GeneDict = _create_gene_dict(genes, chromosome, start, stop)
        self.gene_df = pd.DataFrame.from_dict(self.GeneDict).T

    def add_target_metadata(self, regex_annotations=None, DirectPairDict=False):

        """Added annotations should match the gene target used for downstream analysis."""

        if regex_annotations == None:
            regex_annotations = self.genes
        elif regex_annotations:
            regex_annotations = regex_annotations

        self.guide_df = _add_guide_target_metadata(
            self.guide_df, regex_annotations, DirectPairDict
        )

    def fetch_sequence_context(self, ref_seq_path):

        """"""

        self.GeneDict = _add_context_sequence_to_gene_dict(
            self.GeneDict, self.gene_df, ref_seq_path
        )

    def annotate_position(self):

        """ """
        pass
        #self.guide_df = _annotate_guide_position(self.GeneDict, self.guide_df)
        
        
        
def _annotate_sgRNAs(guides, genes, chrom, start, stop, annotations, DirectPairDict, ref_seq_path):
    
    """Given the metadata csv with targets and Mb positions, annotate the sgRNAs in the screen.genes table."""
    
    
    _guides = _GuideAnnotation(guides, genes, chrom, start, stop)
    _guides.add_target_metadata(annotations, DirectPairDict)
    _guides.fetch_sequence_context(ref_seq_path)
    _guides.annotate_position()
    
    return _guides.guide_df