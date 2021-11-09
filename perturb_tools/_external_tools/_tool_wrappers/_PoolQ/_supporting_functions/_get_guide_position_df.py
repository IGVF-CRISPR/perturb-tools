import re
import pandas as pd

from ..._genome_editing._sequence_manipulation._get_reverse_complement import (
    _get_reverse_complement,
)
from ..._genome_editing._sequence_manipulation._get_chromosome_sequence import (
    _get_chromosome_sequence,
)


def _isolate_target_sequence(
    chromosome,
    start,
    stop,
    ref_seq_path="/home/mvinyard/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
):

    """ """

    # convert from Mb to integer values
    start_ = int(start * 1e6)
    stop_ = int(stop * 1e6)

    target_seq = _get_chromosome_sequence(ref_seq_path, chromosome)[start_:stop_]

    return target_seq


def _find_sgRNA(guide_sequence, target_sequence, strand):

    """
    Find start, stop position of a sgRNA in a target sequence.

    Parameters:
    -----------
    guide_sequence

    target_sequence

    Returns:
    --------
    start, stop
    """

    guide_span = re.search(guide_sequence, target_sequence).span()
    start, stop = guide_span[0], guide_span[1]

    if strand == "-":
        start, stop = (
            len(target_sequence) - guide_span[0],
            len(target_sequence) - guide_span[1],
        )

    return start, stop, strand


def _make_guide_df(guide_sequence, target_sequence):
    """"""

    target_sequence_rc = _get_reverse_complement(target_sequence)

    GuideDict = {}
    nonmatch = []
    #     GuideDict["non-matching"] = []

    for guide in guide_sequence:
        GuideDict[guide] = {}

        try:
            (
                GuideDict[guide]["start"],
                GuideDict[guide]["stop"],
                GuideDict[guide]["strand"],
            ) = _find_sgRNA(guide, target_sequence, strand="+")
        except:
            try:
                (
                    GuideDict[guide]["start"],
                    GuideDict[guide]["stop"],
                    GuideDict[guide]["strand"],
                ) = _find_sgRNA(guide, target_sequence_rc, strand="-")
            except:
                #                 GuideDict["non-matching"].append(guide)
                (
                    GuideDict[guide]["start"],
                    GuideDict[guide]["stop"],
                    GuideDict[guide]["strand"],
                ) = ("nonmatch", "nonmatch", "nonmatch")

    guide_df = (
        pd.DataFrame.from_dict(GuideDict)
        .T.reset_index()
        .rename({"index": "sequence"}, axis=1)
    )

    return guide_df


def _get_guide_position_df(region, chromosome, start, stop, rowsfile, row_delim):

    target_seq = _isolate_target_sequence(chromosome, start, stop)

    guides = pd.read_csv(
        rowsfile, sep=row_delim, header=None, names=["sequence", "target_region"],
    )
    region_guides = guides.loc[guides.target_region.str.contains(region)]
    guide_df = _make_guide_df(
        guide_sequence=region_guides["sequence"], target_sequence=target_seq
    )
    guide_df.start + start
    guide_df.stop + stop

    return guide_df
