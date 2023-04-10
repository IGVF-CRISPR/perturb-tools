import numpy as np


def _log_fold_change(mat: np.ndarray, cond1: int, cond2: int) -> np.ndarray:
    """
    Calculates the log fold-change between two conditions.

    The function takes in an array of log-transformed expression values and indices for two conditions.
    It then calculates the difference in expression between the two conditions (cond1/cond2) for each gene,
    returning an array of log fold-changes.

    Parameters
    ----------
    mat : `numpy.ndarray`
        Array of log-transformed expression values with shape `(n_conditions, n_genes)`.
    cond1 : `int`
        Index of the first condition for calculating the log-fold change.
    cond2 : `int`
        Index of the second condition for calculating the log-fold change.

    Returns
    -------
    log_fc : `numpy.ndarray`
        Array of log fold-changes with shape `(n_genes,)`.
    """
    return mat[cond1, :] - mat[cond2, :]
