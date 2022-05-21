import numpy as np
def _log_fold_change(mat: np.array, cond1, cond2):

    """Log fold-change calculation. Assumes log-transformed values as input."""
    return mat[:, cond1] - mat[:, cond2]
