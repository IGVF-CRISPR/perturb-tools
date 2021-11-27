
# _print_underline.py

__module_name__ = "_print_underline.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import numpy as np

# local imports #
# ------------- #
from ._python_string_formatting import _format_string_printing_font


def _print_underline(string, formatting=["BOLD", "RED"], n_newline=1):

    """
    Format a string and print a matching underline, underneath.
    Parameters:
    -----------
    string [ required ]
        type: str
        
    formatting [ optional ]
        default: ["BOLD", "RED"]
        type: list
        
    n_newline [ optional ]
        default: 1
        type: int
    Returns:
    --------
    None, prints results
    Notes:
    ------
    Examples:
    ---------
    ```
    v.ut.print_underlined("My title:")
    
    >>> My title:
    >>> ---------
    ```
    """

    n_under = len(string)
    underline = []
    for i in range(n_under):
        underline.append("{}".format("-", i))

    formatted_string = _format_string_printing_font("".join(string), formatting)
    formatted_underline = _format_string_printing_font("".join(underline), formatting)

    print(formatted_string)
    print(formatted_underline, "".join(np.repeat("\n", n_newline)))