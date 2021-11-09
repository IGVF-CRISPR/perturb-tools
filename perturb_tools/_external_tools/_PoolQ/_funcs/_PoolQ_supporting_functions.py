from ..._utilities._system_utils._flexible_mkdir import _flexible_multilevel_mkdir
from ..._utilities._ux_utils._pystrings import _format_string_printing_font
from ..._utilities._ux_utils._print_underline import _print_underline
from ..._utilities._data_utils._read_txt import _read_txt

import numpy as np
import pandas as pd
import glob, os
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib

# def _return_outs_filepaths(outsdir, outfile_names):

#     """"""

#     outfiles = []

#     for file in outfile_names:
#         outfiles.append(os.path.join(outsdir, file))

#     return outfiles


def _make_results_dict(poolq_outs_path):

    ResultsDict = {}

    ResultsDict["raw_counts"] = pd.read_csv(
        os.path.join(poolq_outs_path, "counts.txt"), sep="\t"
    )
    ResultsDict["barcode_counts"] = pd.read_csv(
        os.path.join(poolq_outs_path, "barcode-counts.txt"), sep="\t"
    )
    ResultsDict["correlation"] = pd.read_csv(
        os.path.join(poolq_outs_path, "correlation.txt"), sep="\t", index_col=0
    )
    ResultsDict["lognorm_counts"] = pd.read_csv(
        os.path.join(poolq_outs_path, "lognormalized-counts.txt"), sep="\t"
    )

    ResultsDict["poolq3.log"] = _read_txt(os.path.join(poolq_outs_path, "poolq3.log"))
    ResultsDict["quality"] = _read_txt(os.path.join(poolq_outs_path, "quality.txt"))
    ResultsDict["runinfo"] = _read_txt(os.path.join(poolq_outs_path, "runinfo.txt"))
    ResultsDict["unexpected-sequences"] = pd.read_csv(
        os.path.join(poolq_outs_path, "unexpected-sequences.txt"), sep="\t"
    )

    return ResultsDict


def _get_read_files(data_dir, bc):

    """

    Parameters:
    -----------
    data_dir
        directory of poolQ inputs
        type: str

    bc
        string indicator in filename of which fastq contains the barcodes
        type: str

    Returns:
    --------
    col_reads, row_reads

    Notes:
    ------
    (1) list comprehension setup should scale to future, multi-pool analyses.

    """

    fastq_matches = np.array(glob.glob(os.path.join(data_dir, "*.fastq.gz")))
    out = np.array([i for i, v in enumerate(fastq_matches) if bc in v])

    col_reads = barcode = fastq_matches[out][0]  # barcode
    row_reads = read_fq = fastq_matches[~out][0]  # read_fq

    return col_reads, row_reads


def _run_PoolQ(
    data_dir,
    barcode_filename,
    conditions_filename,
    java_path,
    barcode_indicator,
    dry_run,
    run_name,
):

    """
    Execute poolQ commands.

    Parameters:
    -----------
    data_dir
        directory of poolQ inputs
        type: str

    bc
        string indicator in filename of which fastq contains the barcodes
        type: str

    Returns:
    --------
    col_reads, row_reads

    Notes:
    ------
    (1) list comprehension setup should scale to future, multi-pool analyses.

    """

    col_reads, row_reads = _get_read_files(data_dir, barcode_indicator)
    barcodes = os.path.join(data_dir, barcode_filename)
    conditions = os.path.join(data_dir, conditions_filename)

    executable = " ".join(
        [
            "bash",
            java_path,
            "--row-reference {}".format(barcodes),
            "--col-reference {}".format(conditions),
            "--row-reads {}".format(row_reads),
            "--col-reads {}".format(col_reads),
            "--row-barcode-policy PREFIX:CACCG@12 --col-barcode-policy FIXED:0",
        ]
    )

    formatted_runtitle = _format_string_printing_font(
        "Run name: {}".format(run_name), ["BOLD", "RED"]
    )
    _print_underline("Run name: {}".format(run_name))
    _print_underline("PoolQ executable:", formatting=["BOLD", "CYAN"], n_newline=0)
    print(executable)
    print(_format_string_printing_font("-----------------\n", ["BOLD", "CYAN"]))
    if dry_run:
        print(
            _format_string_printing_font(
                "\nDry run complete. Please inspect PoolQ executable.", ["BOLD"]
            )
        )
    else:
        os.system(executable)


def _organize_PoolQ_outputs(run_name, standard_PoolQ_outfiles, outpath, fetch_only):

    """Move the outputs of PoolQ from the default output directory (i.e., the working directory "./") to one defined."""

    run_dir = run_name + ".PoolQ.3.3.2"

    if not fetch_only:
        _flexible_multilevel_mkdir(run_dir)

    for file in standard_PoolQ_outfiles:
        if not fetch_only:
            os.system("mv {} {}".format(file, os.path.join(run_dir)))

    if outpath:
        if not fetch_only:
            os.system("mv {} {}".format(run_dir, outpath))
        outs = os.path.join(outpath, run_dir)
        return outs
    else:
        outs = os.path.join(os.getcwd(), run_dir)
        return outs


def _get_plot_columns(df, ignore_containing):

    """"""

    col_list = []

    for column in df.columns:
        if column in ignore_containing:
            continue
        else:
            col_list.append(column)

    return col_list


def _plot_correlation_scatter(
    df, savename, figsize=1, ignore_containing=["sequence", "barcode_id"]
):

    """
    df
        lognorm_counts_df
    """

    col_list = _get_plot_columns(df, ignore_containing)

    fig = plt.Figure(figsize=(25 * figsize, 25 * figsize))
    gs = GridSpec(len(col_list), len(col_list))
    AxesDict = {}

    for i, col1 in enumerate(col_list):
        AxesDict[i] = {}
        for j, col2 in enumerate(col_list):
            if i + j >= len(col_list):  # CREATES HORIZONTAL
                continue  # CREATES HORIZONTAL
            else:
                AxesDict[i][j] = fig.add_subplot(gs[i, j])
                AxesDict[i][j].scatter(df[col1], df[col2], c="navy", alpha=0.5)
            if i == 0:
                AxesDict[i][j].set_title(col2)
            if j == 0:
                AxesDict[i][j].set_ylabel(col1)
            else:
                continue
    plt.tight_layout()
    fig_savename = savename + ".correlation.scatter.png"
    fig.savefig(fig_savename)
    plt.show()

    return fig, fig_savename
