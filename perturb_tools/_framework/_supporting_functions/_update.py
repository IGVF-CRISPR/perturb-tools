def _update(self, **kwargs):

    for key, value in kwargs.items():
        self.X[key] = value

    return self