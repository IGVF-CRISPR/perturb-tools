def _widen_regions(regions_df, amount):

    regions_df.Start = regions_df.Start - amount
    regions_df.End = regions_df.End + amount

    return regions_df