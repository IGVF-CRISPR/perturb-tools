
def _log_fold_change(df, cond1, cond2):

    """Log fold-change calculation. Assumes log-transformed values as input."""
    return df[cond2] - df[cond1]