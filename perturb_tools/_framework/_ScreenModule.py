
# _ScreenModule.py
__module_name__ = "_ScreenModule.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


import vintools as v


from ._supporting_functions._print_screen_object import _print_screen_object
from ._supporting_functions._data_reading._read_screen_from_PoolQ import _read_screen_from_PoolQ

class _Screen:
    def __init__(self, X=False):

        if X:
            self.X = X
            n_guides, n_conditions = _print_screen_object(self.X)

    def read_PoolQ(self, path):

        """ """

        self._PoolQ_outpath = path
        self._PoolQScreenDict = _read_screen_from_PoolQ(self._PoolQ_outpath)
        
        for key, value in v.ut.update_dict(self._PoolQScreenDict).items():
            self.__setattr__(key, value)
        _print_screen_object(self)