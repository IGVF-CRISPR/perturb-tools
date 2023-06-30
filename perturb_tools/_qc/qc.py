from typing import Optional, List, Union
import numpy as np
import numpy.ma as ma
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData


def _set_sample_correlation(
    screen: AnnData,
    method: str = "pearson",
    count_layer: Optional[str] = None,
    guide_idx: Optional[Union[np.ndarray, List[bool]]] = None,
    prefix: str = "",
):
    """Calculate pairwise correlation between samples in AnnData object and store the result in the object.

    Args:
        screen (AnnData): Annotated data matrix containing gene expression values.
        method (str, optional): Method used for calculating correlation. Defaults to "pearson".
        count_layer (str, optional): Key for the layer to use if AnnData object contains multiple layers. Defaults to None.
        guide_idx (Optional[Union[np.ndarray, List[bool]]], optional): Boolean array or list of indices to subset samples before calculating correlation. Defaults to None.
        prefix (str, optional): Prefix added to name of new columns added to AnnData.obs. Defaults to "".

    Returns:
        None

    Raises:
        ValueError: Raised when an invalid `method` or `count_layer` argument is passed.

    References:
        - https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html
        - https://numpy.org/doc/stable/reference/generated/numpy.ma.corrcoef.html
    """
    corr_func_dict = {
        "pearson": lambda X: ma.corrcoef(X),
        "spearman": lambda X: scipy.stats.spearmanr(X.T, nan_policy="omit")[0],
    }
    corr_func = corr_func_dict[method]
    if guide_idx is not None:
        screen_subset = screen[:, guide_idx].copy()
    else:
        screen_subset = screen
    if count_layer is None:
        count_matrix = screen_subset.X
        count_layer_label = "X"
    else:
        if count_layer not in screen_subset.layers:
            raise ValueError(f"Invalid layer key {count_layer}.")
        count_matrix = screen_subset.layers[count_layer]
        count_layer_label = count_layer
    if prefix != "":
        prefix = f"{prefix}_"
    screen.obsm[f"{prefix}corr_{count_layer_label}"] = corr_func(
        ma.masked_invalid(count_matrix)
    )

    screen.obs[f"{prefix}mean_corr_{count_layer_label}"] = screen.obsm[
        f"{prefix}corr_{count_layer_label}"
    ].mean(0)

    screen.obs[f"{prefix}median_corr_{count_layer_label}"] = np.median(
        screen.obsm[f"{prefix}corr_{count_layer_label}"], axis=1
    )


def get_outlier_guides(
    screen: AnnData,
    cond_col: str,
    count_layer: Optional[str] = None,
    mad_z_thres: float = 5,
    abs_RPM_thres: float = 10000,
) -> pd.DataFrame:
    """Obtain outlier guides.

    Args:
        screen (AnnData): Annotated data matrix.
        cond_col (str): Column name in obs where experimental conditions are annotated.
        count_layer (Optional[str], optional): Name of layer used for quantification. Defaults to None.
        mad_z_thres (float, optional): Threshold (in MAD units) above which guides will be considered as outliers. Defaults to 5.
        abs_RPM_thres (float, optional): RPM threshold above which guides will be considered as outliers. Defaults to 10000.

    Returns:
        pd.DataFrame: A pandas DataFrame with a list of outlier guides per experimental condition.

    Raises:
        ValueError: If `count_layer` is not found in `screen.layers`

    This function obtains outlier guides by comparing the counts of guides within experimental conditions to each other.
    Outlier guides are defined as those that show extreme counts compared to guides based on a specified threshold (`mad_z_thres`). Guides should also have a high absolute RPM value (`abs_RPM_thres`) if they are to be considered as outliers.
    """
    if count_layer is None:
        count_layer_label = "X"
        count_matrix = screen.X
    else:
        if count_layer not in screen.layers:
            raise ValueError(f"{count_layer} not in screen.layers")
        count_layer_label = count_layer
        count_matrix = screen.layers[count_layer]
    if f"{count_layer_label}_RPM" not in screen.layers:
        screen.layers[f"{count_layer_label}_RPM"] = (
            count_matrix / count_matrix.sum(axis=1)[:, None] * 1e6
        )
    aberr_dict = {}

    for cond in screen.obs[cond_col].unique():
        screen_subset = screen[screen.obs[cond_col] == cond, :].copy()
        median_p = np.nanmedian(
            screen_subset.layers[f"{count_layer_label}_RPM"], axis=1
        )
        aberr_guide_df_condit = []
        for i, sample in enumerate(screen_subset.obs.index):
            outlier_idx = np.where(
                (
                    screen_subset.layers[f"{count_layer_label}_RPM"][:, i]
                    > median_p * mad_z_thres
                )
                & (
                    screen_subset.layers[f"{count_layer_label}_RPM"][:, i]
                    > abs_RPM_thres
                )
            )[0]
            outlier_guides = screen_subset.var.iloc[outlier_idx, :].copy()
            outlier_guides["sample"] = screen_subset.obs.index[i]
            outlier_guides["RPM"] = screen_subset.layers[f"{count_layer_label}_RPM"][
                outlier_idx, i
            ]
            aberr_guide_df_condit.append(outlier_guides[["sample", "RPM"]])
        aberr_dict[cond] = aberr_guide_df_condit
    aberr_guide_dfs = []
    for df in aberr_dict.values():
        aberr_guide_dfs.extend(df)
    aberr_guides = pd.concat(aberr_guide_dfs, axis=0)
    # aberr_guides.index = aberr_idx_list
    return aberr_guides.reset_index()
