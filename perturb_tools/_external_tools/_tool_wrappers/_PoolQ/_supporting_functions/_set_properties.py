class _PropertiesDict(object):
    def __init__(self, PathDict):

        self.PathDict = PathDict
        self.verbose = True
        self.java_path = os.path.join(
            v.__path__[0], "_tools/_PoolQ/poolq-3.3.2/poolq3.sh"
        )
        self.barcode_indicator = "barcode"
        self.barcode_filename = "rows.txt"
        self.conditions_filename = "column.txt"
        self.standard_PoolQ_outfiles = [
            "barcode-counts.txt",
            "correlation.txt",
            "counts.txt",
            "lognormalized-counts.txt",
            "poolq3.log",
            "quality.txt",
            "runinfo.txt",
            "unexpected-sequences.txt",
        ]

    def _update_parameter(self, parameter, value):

        self.__setattr__(parameter, value)

    def update(self, **kwargs):

        """
        Update all passed parameters in a ParamDict.
        """

        for parameter, value in kwargs.items():
            self._update_parameter(parameter, value)
            message1 = v.ut.format_pystring("Parameter adjusted: ", ["BOLD", "RED"])
            message2 = v.ut.format_pystring(
                "{} = {}\n".format(parameter, value), ["BOLD"]
            )
            print(message1, message2)


def _set_properties(PathDict, **kwargs):

    """"""

    PropertiesDict = _PropertiesDict(PathDict)
    PropertiesDict.update(**kwargs)

    return PropertiesDict