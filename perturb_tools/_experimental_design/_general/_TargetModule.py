
from ._GTF_Module import _GTF
from ._get_screen_target_features import _gene_from_gtf, _get_feature, _get_gene_body_bounds
from ._identify_overlapping_fragments import _OverlappingFragments
from ._query_region_for_PAM import _query_region_for_PAM

from ._make_guide_library_df import _make_guide_library_df
from ._annotate_protospacer import _annotate_protospacer
from ._add_sgRNA_context import _add_sgRNA_context
from ._scan_for_BsmbI import _scan_for_BsmbI

from ._widen_regions import _widen_regions
from ._add_region_labels import _add_region_labels

class _ScreeningTarget:
    
    def __init__(self, GTF, ref_seq_path):
        
        self.ref_seq_path = ref_seq_path
        
        if type(GTF) == str:
            gtf = _GTF(GTF)
            self.gtf_file = gtf.file
        else:
            self.gtf_file = GTF.file
        
    def get_gene(self, chromosome, gene):
        
        self.gene = gene
        self.gene_df = _gene_from_gtf(self.gtf_file, chromosome, gene)
        self.gene_min, self.gene_max = _get_gene_body_bounds(self.gene_df)
        
    def get_features(self, feature):
        
        self.feature_df = _get_feature(self.gene_df, feature)

    def get_exons(self):
        
        
        self.exon_df = _get_feature(self.gene_df, feature="exon")
        self.exon_df = self.exon_df[[0, 3, 4]]
        self.exon_df.columns = ["Chromosome", "Start", "End"]
        
#         print("\nCheck the exon sequences for overlaps...\n")
        _exon_overlaps = _OverlappingFragments(self.exon_df)
        _exon_overlaps.find_all()
        
        self.exon_df = _exon_overlaps.df
        self.exon_df['Start'] = self.exon_df['Start'].astype(int)
        self.exon_df['End'] = self.exon_df['End'].astype(int)
        self.exon_df['exon_length'] = self.exon_df["End"].astype(int) - self.exon_df["Start"].astype(int)
        
        
    def build_guide_library(self, exons=True, feature=False, PAM="NGG", filter_BsmbI=True, widen=False):
        
        """
        forward_pams, reverse_pams, KMT2C, chrom_seqs
        """
        
        if exons == True:
            print("\nUsing {} exons...\n".format(self.gene))
            self.regions_df = self.exon_df
        else:
            try:
                self.regions_df = self.feature_df
            except:
                self.regions_df = self._get_feature(self.gene_df, feature)
        if widen:
            print("\nWidening regions by {} nucleotides on either side of each boundary...\n".format(widen))
            self.regions_df = _widen_regions(self.regions_df, widen)
                
        print("\nBuilding sgRNA library DataFrame...")
        
        forward_pams, reverse_pams, chrom_seq = _query_region_for_PAM(self.regions_df, self.ref_seq_path, PAM=PAM)
        
        print("PAMs, Chromosome Seqeunce isolated")
        
        self.guide_df = _make_guide_library_df(forward_pams, reverse_pams, self.regions_df, chrom_seq)
        
        print("\nAdding sgRNA sequence context")
        self.guide_df = _add_sgRNA_context(self.guide_df, chrom_seq)
#         self.guide_df = _add_region_labels(self.guide_df)
                
        print("\nAdding BsmbI annotations")
        self.guide_df, self.BsmbI_guide_df = _scan_for_BsmbI(self.guide_df, filter_BsmbI)
        self.guide_df = self.guide_df.reset_index(drop=True)
        
    def write(filename=False):
        
        if not filename:
            filename = "perturb_tools.target_library.{}.xlsx".format(self.gene)

        self.guide_df.to_excel(filename)