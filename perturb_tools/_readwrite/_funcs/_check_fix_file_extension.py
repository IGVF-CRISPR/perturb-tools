# _check_fix_file_extension.py
__module_name__ = "_check_fix_file_extension.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


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