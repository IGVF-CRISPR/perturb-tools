from typing import Union, List, Literal
import scipy.sparse
import numpy as np
import mudata as md


def _check_args(mdata, mod1, mod2, *args):
    for arg in args:
        if not arg in mdata[mod1].var.columns:
            raise ValueError(
                f"{arg} not in mdata[{mod1}].var.columns. Check the input."
            )
        if not arg in mdata[mod2].var.columns:
            raise ValueError(
                f"{arg} not in mdata[{mod2}].var.columns. Check the input."
            )


def calculate_same_chrom(
    mdata: md.MuData,
    mod1: str,
    mod2: str,
    output_key: str = "same_chrom",
    chrom_col: Union[List[str], str] = "chr",
):
    _check_args(
        mdata,
        mod1,
        mod2,
    )
    same_chrom = (
        mdata[mod1].var[chrom_col].astype(str).values[:, None]
        == mdata[mod2].var[chrom_col].astype(str)[None, :]
    )
    assert same_chrom.shape == (mdata[mod1].n_var, mdata[mod2].n_var)
    full_varp_mat = np.zeros((mdata.n_vars, mdata.n_vars))
    mod1_idx = mdata.varmap[mod1]
    mod1_idx = mod1_idx[mod1_idx > 0] - 1
    mod2_idx = mdata.varmap[mod2]
    mod2_idx = mod2_idx[mod2_idx > 0] - 1
    full_varp_mat[
        np.repeat(mod1_idx, len(mod2_idx)), np.tile(mod2_idx, len(mod1_idx))
    ] = same_chrom.flat
    mdata.varp[output_key] = scipy.sparse.coo_matrix(full_varp_mat)


def calculate_distance(
    mdata: md.MuData,
    mod1: str,
    mod2: str,
    same_chrom_key: str = "same_chrom",
    output_key: str = "distance",
    chrom_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
    how: Literal["midpoint", "shortest"] = "midpoint",
):
    """Calculate the pairwise distance between .var entries in mod1 and mod2
    and saves the result in mdata.varp[output_key]
    """
    _check_args(mdata, mod1, mod2, chrom_col, start_col, end_col)
    if same_chrom_key not in mdata.varp.keys():
        calculate_same_chrom(mdata, mod1, mod2, same_chrom_key, chrom_col)
    distance = np.abs(mdata.var[f"{mod1}:{start_col}"].astype(float).values[:,None] - mdata.var[f"{mod2}:{start_col}"].astype(float).values)
    distance = distance * mdata.varp[same_chrom_key]
    mdata.varp[output_key] = scipy.sparse.coo_matrix(distance)
    
