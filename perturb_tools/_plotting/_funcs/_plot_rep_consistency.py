import seaborn as sns
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

def _annotate_pearsonr_foldchange(x, y, **kws):
    ignore_idx = np.isnan(x) | np.isnan(y)
    x = x[~ignore_idx]
    y = y[~ignore_idx]
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    
    # count how many annotations are already present
    n = len([c for c in ax.get_children() if isinstance(c, matplotlib.text.Annotation)])
    pos = (.1, .9 - .1*n)
    
    ax.annotate("r={:.2f}".format(r),
                xy=pos, xycoords=ax.transAxes)

def plot_rep_consistency_lfc(
        screen,
        cond1,
        cond2,
        hue = None,
        lognorm_counts_key = "lognorm_counts",
        rep_condit = "replicate",
        compare_condit = "sort",
        out_guides_suffix = "lfc",
        keep_rep_lfc = False,
        title = "Replicate consistency of LFC",
        return_plot = True,
        **kws
        ):

    if not hue is None and not hue in screen.guides.columns:
        raise KeyError("{} not in guide features.".format(hue))

    lfcs_df_columns =  [s + ".{}_{}.{}".format(cond1, cond2, out_guides_suffix) for s in screen.condit[rep_condit].unique()]
    lfcs_df = screen.log_fold_change_reps(
        cond1,
        cond2,
        lognorm_counts_key = lognorm_counts_key,
        rep_condit = rep_condit,
        compare_condit = compare_condit,
        out_guides_suffix = out_guides_suffix,
        keep_result = keep_rep_lfc,
        )

    if not hue is None:
        lfcs_df = pd.concat((screen.guides[[hue]], lfcs_df), axis = 1)    
    try:
        pairplot_obj = sns.pairplot(lfcs_df, kind = "scatter", hue = hue, corner = True,
                **kws)
    except:
        return(lfcs_df)
    pairplot_obj.map_lower(_annotate_pearsonr_foldchange)
    #pairplot_obj.fig.suptitle(title)
    #pairplot_obj.fig.tight_layout()
    if return_plot:
        return(pairplot_obj)
