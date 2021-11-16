
import os, glob
import numpy as np

from ...._utilities._funcs._python_string_formatting import _format_string_printing_font
from ...._utilities._funcs._print_underline import _print_underline

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
    
    print("Datadir:", data_dir)

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