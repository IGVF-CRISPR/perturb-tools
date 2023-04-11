# import plotly.express as px
from tqdm import tqdm_notebook as tqdm

#from ..._utilities._funcs._get_gene_exons import _get_gene_exons

from ._plot_gene_exon import _plot_gene_exon

def _plot_guide_outline(fig, start, end, y):

    """"""

    fig.add_shape(
        x0=start,
        x1=end,
        y0=0,
        y1=y,
        line_width=2,
    )


    
def _plot_guide_fill(fig, start, end, y, fillcolor="navy", opacity=0.1):

    """"""

    fig.add_shape(
        x0=start,
        x1=end,
        y0=0,
        y1=y,
        line_width=0,
        fillcolor=fillcolor,
        opacity=opacity,
    )



def _plot_single_guide(fig, start, end, y, fillcolor="navy", opacity=0.1):

    _plot_guide_outline(fig, start, end, y)
    _plot_guide_fill(fig, start, end, y, fillcolor="navy", opacity=0.1)


def _format_plot(fig, 
                 experiment_guide_df, 
                 y_title="\u0394Log\u2082(Fold Change)", 
                 gridcolor="mistyrose"):

    """"""
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0, 0, 0, 0)",
            "paper_bgcolor": "rgba(0, 0, 0, 0)",
        }
    )
    fig.add_hline(y=0.0, line_width=0.5, line_dash="dash", line_color="red")
    fig.update_layout(font_family="Avenir", title_font_family="sans-serif")
    fig.update_yaxes(title=y_title)
    fig.update_xaxes(
        title="Genomic Coordinates: {}".format(experiment_guide_df.Chromosome.unique()[0])
    )
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor=gridcolor)
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor=gridcolor)


def _plot_guide_enrichment(experiment_guide_df, title=None, return_fig=False, exon_color='red', gene=False, gtf_path=None):

    """

    Minimal func for now. Lacks several helpful / usability things and is very stiff.

    """

    fig = px.scatter(
        experiment_guide_df,
        title=title,
        x="center",
        y="lfc.mean",
    )

    _format_plot(fig, experiment_guide_df)

    for i, row in tqdm(experiment_guide_df.iterrows()):
        _plot_single_guide(fig, row["Start"], row["End"], row["lfc.mean"])
    
    if gene:
        if type(gene) == str:
            gene = [gene]
        for i in range(len(gene)):
            gene_df = _get_gene_exons(gene[i], experiment_guide_df.Chromosome.unique()[0], gtf_path)
            _plot_gene_exon(fig, gene[i], gene_df, lfc=experiment_guide_df['lfc.mean'], exon_color=exon_color)
    fig.show()

    if return_fig:
        return fig
