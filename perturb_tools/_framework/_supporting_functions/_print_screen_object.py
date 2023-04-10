# _print_screen_object.py
__module_name__ = "_print_screen_object.py"
__author__ = ", ".join(["Michael E. Vinyard", "Jayoung Ryu"])
__email__ = ", ".join(["vinyard@g.harvard.edu", "jayoung_ryu@g.harvard.edu"])
from anndata import AnnData


def _print_screen_object(ScreenObject: AnnData):
    """
    Would be good to return an organizational dictionary
    """

    n_samples, n_guides = ScreenObject.shape[0], ScreenObject.shape[1]

    descr = (
        "Genome Editing Screen comprised of n_samples x n_guides = {} x {}\n".format(
            n_guides, n_samples
        )
    )

    MainScreenAttributes = [
        "guides",
        "condit",
        "condit_m",
        "condit_p",
        "layers",
        "uns",
    ]

    for attribute in MainScreenAttributes:
        descr += "   {: <11}".format(attribute + ":")
        x = list(ScreenObject.__getattribute__(attribute).keys())
        for n, i in enumerate(x):
            if n != len(x) - 1:
                descr += "'{}', ".format(i)
            else:
                descr += "'{}'".format(i)
        descr += "\n"

    return n_guides, n_samples, descr
