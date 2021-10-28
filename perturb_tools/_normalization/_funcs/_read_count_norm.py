
import numpy as np

def _read_count_normalize(X):

    """Read depth normalization by sample. Assumes samples are columns and guides are rows."""

    return (X / X.sum(axis=0)) * 1e6


def _log_transform_read_count(X):

    """"""
    return np.log2(X + 1)


def _log_normalize_read_count(X):
    
    """Following the protocol written clearly, here: 
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0170445#sec002 
    (see Methods). 
    """
    
    X_read_norm = _read_count_normalize(X)
    X_log_read_norm = _log_transform_read_count(X_read_norm)

    return X_log_read_norm