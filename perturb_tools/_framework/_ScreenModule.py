# _ScreenModule.py
__module_name__ = "_ScreenModule.py"
__author__ = ", ".join(["Michael E. Vinyard", "Jayoung Kim Ryu"])
__email__ = ", ".join(["vinyard@g.harvard.edu", "jayoung_ryu@g.harvard.edu"])

from anndata import AnnData

# from ._supporting_functions._guides._GuideAnnotationModule import _annotate_sgRNAs


def annotate_guides(
    adata: AnnData, genes, chrom, start, stop, annotations, DirectPairDict, ref_seq_path
):
    """
    Annotate sgRNA table.

    """
    # adata.var = _annotate_sgRNAs(
    #     adata.var,
    #     genes,
    #     chrom,
    #     start,
    #     stop,
    #     annotations,
    #     DirectPairDict,
    #     ref_seq_path,
    # )

    NotImplemented
