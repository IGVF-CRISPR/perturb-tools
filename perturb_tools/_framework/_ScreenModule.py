# _ScreenModule.py


class Screen:
    def __init__(self, X=False):

        if X:
            self.X = X
            n_guides, n_conditions = _report_func(self.X)

    def read_PoolQ(self, path):

        """ """

        self._PoolQ_outpath = path
        self._PoolQScreenDict = _read_screen_from_PoolQ(self._PoolQ_outpath)

        for key, value in _update(self._PoolQScreenDict).items():
            self.__setattr__(key, value)

        _report_func(self)