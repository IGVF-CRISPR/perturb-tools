
import numpy as np

def _determine_overlap_3prime(
    fragment1_start, fragment1_stop, fragment2_start, fragment2_stop
):

    """
    Observe whether the start-site of fragment2 overlaps the 3' end of fragment1.

    Parameters:
    -----------
    fragment1_start

    fragment1_stop

    fragment2_start

    fragment2_stop

    Returns:
    --------
    overlap_3prime
        type: bool

    Notes:
    ------
    """

    return (fragment2_start < fragment1_stop) and (fragment2_stop > fragment1_stop)


def _determine_overlap_5prime(
    fragment1_start, fragment1_stop, fragment2_start, fragment2_stop
):

    """
    Observe whether the start-site of fragment2 overlaps the 5' end of fragment1.

    Parameters:
    -----------
    fragment1_start

    fragment1_stop

    fragment2_start

    fragment2_stop

    Returns:
    --------
    overlap_5prime
        type: bool

    Notes:
    ------
    """

    return (fragment2_start < fragment1_start) and (fragment2_stop > fragment1_start)


def _determine_partial_overlapping_fragments(df, silent, start_key="Start", end_key="End"):

    """

    Parameters:
    -----------
    df
        requires the following columns: ['start', 'stop']

    Returns:
    --------
    OverlapDict
    """

    OverlapDict = {}
    OverlapDict["5'"] = []
    OverlapDict["3'"] = []

    for fragment1 in range(len(df)):
        for fragment2 in range(len(df)):

            if _determine_overlap_3prime(
                df[start_key].iloc[fragment1],
                df[end_key].iloc[fragment1],
                df[start_key].iloc[fragment2],
                df[end_key].iloc[fragment2],
            ):
                OverlapDict["5'"].append([fragment1, fragment2])
                if not silent:
                    print(
                        "fragment {} overlaps the 5' end of fragment {}".format(
                            fragment1, fragment2
                        )
                    )

            if _determine_overlap_5prime(
                df[start_key].iloc[fragment1],
                df[end_key].iloc[fragment1],
                df[start_key].iloc[fragment2],
                df[end_key].iloc[fragment2],
            ):
                OverlapDict["3'"].append([fragment1, fragment2])
                if not silent:
                    print(
                        "fragment {} overlaps the 3' end of fragment {}".format(
                            fragment1, fragment2
                        )
                    )
    return OverlapDict


def _determine_total_overlap(
    fragment1_start, fragment1_stop, fragment2_start, fragment2_stop
):
    return (fragment1_start < fragment2_start) and (fragment1_stop > fragment2_stop)


def _determine_total_overlapping_fragments(df, silent=False, start_key="Start", end_key="End"):

    overlapped_fragments = []

    for fragment1 in range(len(df)):
        for fragment2 in range(len(df)):

            if _determine_total_overlap(
                df[start_key].iloc[fragment1],
                df[end_key].iloc[fragment1],
                df[start_key].iloc[fragment2],
                df[end_key].iloc[fragment2],
            ):
                overlapped_fragments.append(fragment1)
                if not silent:
                    print("fragment {} engulfs {}".format(fragment1, fragment2))
            elif _determine_total_overlap(
                df[start_key].iloc[fragment2],
                df[end_key].iloc[fragment2],
                df[start_key].iloc[fragment1],
                df[end_key].iloc[fragment1],
            ):
                overlapped_fragments.append(fragment2)
                if not silent:
                    print("fragment {} engulfs {}".format(fragment2, fragment1))

    return np.unique(overlapped_fragments)


def _drop_totally_overlapped_fragments(df, overlapped_fragments_list, start_key="Start", end_key="End"):

    """"""

    print("\nDropping the following exonic regions due to overlap:\n")
    for i in overlapped_fragments_list:
        print("\t{}-{}".format(df[start_key].iloc[i], df[end_key].iloc[i]))

    df_ = df.drop(overlapped_fragments_list)

    return df_.reset_index(drop=True)


def _identify_unique_fragments(df, start_key="Start", end_key="End"):

    df = df.reset_index(drop=True)
    df_ = (
        df.sort_values(start_key)
        .reset_index(drop=True)
        .drop_duplicates(subset=end_key, keep="first")
    )

    _df_ = (
        df_.sort_values(end_key)
        .reset_index(drop=True)
        .drop_duplicates(subset=start_key, keep="first")
    )

    return _df_.reset_index(drop=True)

class _OverlappingFragments:
    def __init__(self, df, silent=False):

        self.silent = silent
        if not self.silent:
            print("Removing duplicate rows")
        self.df = df.drop_duplicates()
        self.df = _identify_unique_fragments(df)
        self.df_raw = df
        self.OverlappingFragmentDict = {}

    def find_total(self):

        self.OverlappingFragmentDict["total"] = _determine_total_overlapping_fragments(
            self.df, self.silent
        )

        self.df = _drop_totally_overlapped_fragments(
            self.df, self.OverlappingFragmentDict["total"]
        )

    def find_equal_endpoint_fragments(self):

        self.OverlappingFragmentDict["longer_idx"] = _identify_unique_short_fragments(
            self.df
        )

        self.df = self.df.iloc[self.OverlappingFragmentDict["longer_idx"]]

    def find_partial(self):

        self.OverlappingFragmentDict[
            "partial"
        ] = _determine_partial_overlapping_fragments(self.df, self.silent)

    def find_all(self):

        self.OverlappingFragmentDict["total"] = _determine_total_overlapping_fragments(
            self.df, self.silent
        )

        self.df = _drop_totally_overlapped_fragments(
            self.df, self.OverlappingFragmentDict["total"]
        )
        self.OverlappingFragmentDict[
            "partial"
        ] = _determine_partial_overlapping_fragments(self.df, self.silent)