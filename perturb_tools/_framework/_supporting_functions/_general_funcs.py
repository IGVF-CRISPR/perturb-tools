def _get_basename_without_extension(path):

    try:
        return os.path.basename(path).split(".")[0]
    except:
        return os.path.basename(path)


def _check_fix_file_extension(filepath, extension, silent=False):

    if not filepath.endswith(extension):
        if not silent:
            print(
                "filepath: {} does not have the correct extension: {}\nAppending: {}...".format(
                    filepath, extension, extension
                )
            )
        return filepath + extension
    else:
        return filepath
    
    
def _create_EmptyDict(keys):

    """Create an empty python Dictionary given a set of keys."""

    Dict = {}
    for key in keys:
        Dict[key] = {}

    return Dict


def _fix_path_for_glob(path):

    if path.endswith("/*"):
        pass
    elif path.endswith("/"):
        path = path + "*"
    else:
        path = path + "/*"

    return path


def _glob_dict(path, verbose=False):

    """
    Loads filepaths from a glob'd directory.

    Parameters:
    -----------
    path
        path to directory containing files.
        type: str

    verbose
        If true, prints FileHandler messages.
        default: False
        type: bool

    Returns:
    --------
    Dict
        keys: filenames
        values: files loaded into memory

    Notes:
    ------
    (1) Can pass with or without "/*"
    """

    filepaths = glob.glob(_fix_path_for_glob(path))

    Dict = {}
    Dict["paths"] = []

    for path in filepaths:
        Dict["paths"].append(path)
        name = _get_basename_without_extension(path)
        file = FileHandler(filepath=path, verbose=verbose)
        Dict[name] = file.read(return_file=True)

    return Dict


def _write_df_to_excel(
    df_list, workbook_path, sheetnames=False, index=False, silent=False
):

    workbook_path = _check_fix_file_extension(workbook_path, ".xlsx", silent)

    if sheetnames:
        assert len(df_list) == len(sheetnames), print(
            "Length of df_list must equal the length of sheet names passed."
        )
    with pd.ExcelWriter(workbook_path) as writer:

        if not silent:
            print("Writing to: {}\n".format(workbook_path))

        for n, df_ in enumerate(df_list):
            if not sheetnames:
                sheet = "Sheet_{}".format(n)
            else:
                sheet = sheetnames[n]
            if not silent:
                print("\tSheet {}:\t{}".format(str(int(n + 1)), sheet))
            df_.to_excel(writer, sheet, index=index)