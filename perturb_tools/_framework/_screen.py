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


def add(screen1: AnnData, screen2: AnnData) -> AnnData:
    """Add two AnnData matrices if their .var and .obs index matches.

    Args:
        screen1 (AnnData): First matrix to add.
        screen2 (AnnData): Second matrix to add.

    Returns:
        AnnData: Resultant matrix after adding screen1 and screen2.

    Raises:
        ValueError: If the row and column indices of both matrices don't match.
    """
    if all(screen1.var.index == screen2.var.index) and all(
        screen1.obs.index == screen2.obs.index
    ):
        return AnnData(
            X=screen1.X + screen2.X, var=screen1.var.copy(), obs=screen1.obs.copy()
        )
    else:
        raise ValueError("Guides/sample description mismatch")
