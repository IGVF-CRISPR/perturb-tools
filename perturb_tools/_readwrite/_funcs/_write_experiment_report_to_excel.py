# _write_experiment_report_to_excel.py

__module_name__ = "_write_experiment_report_to_excel.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import pandas as pd


# local imports #
# ------------- #
from ._check_fix_file_extension import _check_fix_file_extension
from ._write_screen_to_excel import _collect_screen_dfs


def _write_experiment_report_to_excel(
    screen,
    experiment_df_list,
    experiment_sheet_names,
    workbook_path="CRISPR.screen.experiment.summary.xlsx",
    silent=False,
    index=False,
):

    df_list, sheetnames = _collect_screen_dfs(screen)
    workbook_path = _check_fix_file_extension(workbook_path, ".xlsx", silent)

    for i in range(len(experiment_df_list)):
        df_list.append(experiment_df_list[i])
        sheetnames.append(experiment_sheet_names[i])

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