# _update_dict.py

__module_name__ = "_update_dict.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


def _update_dict(self, **kwargs):
    
    """Update a dictionary given some keyword argument pair."""

    for key, value in kwargs.items():
        self.X[key] = value

    return self
