
def _fold_change(df, cond1, cond2):

    """Linear fold-change calculation. Assumes non-log-transformed values as input."""
    return df[cond2] / df[cond1]