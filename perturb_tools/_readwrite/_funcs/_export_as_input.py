from typing import Optional, Union
import pandas as pd

from anndata import AnnData


def to_mageck_input(
    adata: AnnData,
    out_path: Optional[str] = None,
    count_layer: Optional[str] = None,
    sgrna_column: Optional[str] = None,
    target_column: str = "target_id",
    sample_prefix: str = "",
) -> Union[pd.DataFrame, None]:
    """
    Formats AnnData object into MAGeCK input format.

    This function converts the gene expression data in AnnData object into a matrix with sgRNAs as
    rows and samples as columns, where each entry represents the count of a unique sgRNA-sequence in a sample.
    The resulting matrix is then saved as a tab-separated file in MAGECK input format or returned as a DataFrame.

    Parameters
    ----------
    adata : `anndata.AnnData`
        Annotated data matrix containing the gene expression data.
    out_path : `str`, optional (default: None)
        Path to save the formatted output file in MAGECK input format. If not provided, returns a dataframe.
    count_layer : `str`, optional (default: None)
        Key to access count data in `adata`. If not specified, uses the default count layer `adata.X`.
    sgrna_column : `str`, optional (default: None)
        Column name from within `adata.var` that corresponds to sgRNA sequences. If not specified, assumes sgRNAs
        are stored as index names in `adata.var` DataFrame.
    target_column : `str`, optional (default: "target_id")
        Column name from within `adata.var` that corresponds to gene names.
    sample_prefix : `str`, optional (default: "")
        Prefix string to insert before each sample name in the resulting matrix.

    Returns
    -------
    mageck_input_df : `pandas.DataFrame` or None
        DataFrame containing MAGECK input-formatted data or None if `out_path` provided.
    """
    if count_layer is None:
        count_matrix = adata.X.T
    else:
        try:
            count_matrix = adata.layers[count_layer]
        except KeyError as exc:
            raise KeyError(
                f"Layer {count_layer} doesn't exist in AnnData object with layers {adata.layers.keys()}"
            ) from exc
    mageck_input_df = (
        pd.DataFrame(
            count_matrix,
            columns=sample_prefix + adata.obs.index,
            index=adata.var.index,
        )
        .fillna(0)
        .astype(int)
    )
    if sgrna_column is None:
        mageck_input_df.insert(0, "sgRNA", adata.var.index.tolist())
    elif sgrna_column in adata.var.columns:
        mageck_input_df.insert(0, "sgRNA", adata.var[sgrna_column])
    elif adata.var.index.name == sgrna_column:
        mageck_input_df.insert(0, "sgRNA", adata.var.index.tolist())
    else:
        raise ValueError(f"{sgrna_column} not found in AnnData.guides.")
    mageck_input_df["sgRNA"] = mageck_input_df["sgRNA"].map(
        lambda s: s.replace(" ", "_")
    )
    mageck_input_df.insert(1, "gene", adata.var[target_column])
    mageck_input_df = mageck_input_df.loc[
        (mageck_input_df.gene.map(lambda o: not pd.isnull(o)))
        & mageck_input_df.gene.map(bool),
        :,
    ]
    if out_path is None:
        return mageck_input_df
    else:
        mageck_input_df.to_csv(out_path, sep="\t", index=False)
