__module_name__ = "_read_screen_from_csv.py"
import pandas as pd
from perturb_tools import Screen
def _read_csv(X_path=None,guide_path=None,condit_path=None,sep=","):
  if not X_path is None:
    X_df = pd.read_csv(X_path, delimiter=sep, header=0, index_col=0)
    X = X_df.values
  else: X=None
  if not guide_path is None:
    guide_df = pd.read_csv(guide_path, sep=sep)
  else: guide_df=None
  if not condit_path is None:
    condit_df = pd.read_csv(condit_path, sep=sep)
  else: condit_df=None

  return Screen(X=X, guides=guide_df, condit=condit_df)