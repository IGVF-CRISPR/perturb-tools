# _write_multi_df_to_excel.py
__module_name__ = "_write_multi_df_to_excel.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# package imports #
# --------------- #
import pandas as pd

# local imports #
# ------------- #
from ._check_fix_file_extension import _check_fix_file_extension


def _append_dfs_and_name_from_dict(df_dict, df_list, sheet_names):
    
    subclass = "df_dict={}".format(df_dict).split("=")[0]
    
    for key in df_dict.keys():
        df_list.append(df_dict[key])
        sheet_names.append(".".join([subclass, key]))
    return df_list, sheet_names

def _collect_screen_dfs(screen):

    """collect all main elements of a screen as a pandas DataFrame."""

    df_list = [pd.DataFrame(screen.X, index = screen.guides.index, columns = screen.condit.index)]
    sheet_names = ["X"]
    for key, mat in screen.layers.items():
        df_list.append(pd.DataFrame(mat, index = screen.guides.index, columns = screen.condit.index))
        sheet_names.append(key)

    df_list.extend([screen.guides, screen.condit])
    sheet_names.extend(['guides', 'condit'])

    for i in [screen.condit_m, screen.condit_p]:
        df_list, sheet_names = _append_dfs_and_name_from_dict(i, df_list, sheet_names)

    
    for key in screen.uns.keys():
        df_list.append(pd.DataFrame(screen.uns[key]))
        sheet_names.append(".".join(["screen.uns", key]))
        
    return df_list, sheet_names


def _write_screen_to_excel(
    screen, workbook_path="./CRISPR_screen.xlsx", index=True, silent=False
):
    
    """
    Write one or more pandas.DataFrames to individual sheets of an Excel Workbook (workbook.xlsx). 
    
    Parameters:
    -----------
    df_list
        type: list(pd.DataFrame, ...)
    
    workbook_path
        type: str
   
    index
        default: False
        type: bool
        
    silent
        default: False
        type: bool
        
    Returns:
    --------
    None
        writes to excel workbook (.xlsx)
        
    Notes:
    ------
    (1) 
    """
    
    df_list, sheetnames = _collect_screen_dfs(screen)
    workbook_path = _check_fix_file_extension(workbook_path, ".xlsx", silent)
    
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
