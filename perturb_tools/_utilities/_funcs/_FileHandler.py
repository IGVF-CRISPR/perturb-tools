
# _FileHandler.py
__module_name__ = "_FileHandler.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

class _FileHandler:
    
    def __init__(self, filepath=False, verbose=True):
        
        """Uses only pure python"""
        
        self.verbose = verbose
        
        if filepath:
            self.filepath = filepath
    
    def read(self, return_file=False, filepath=False):
        
        """
        """
        
        if filepath:
            self.filepath = filepath
        if self.verbose:
            print("Reading: {}".format(self.filepath))
        with open(self.filepath, "r") as f:
            self.file = f.readlines()
        
        if return_file:
            return self.file
        
    def count(self):
        
        self.file_length = len(self.file)
        if self.verbose:
            print(self.file_length)