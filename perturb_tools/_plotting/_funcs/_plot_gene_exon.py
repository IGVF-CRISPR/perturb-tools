
import numpy as np

def _get_gene_y_position(lfc):
    return lfc.min() - (lfc.max() - lfc.min()) / 10

def _plot_gene_single_exon(fig, start, end, y, fillcolor="red", opacity=0.8):

    """"""

    fig.add_shape(
        x0=start,
        x1=end,
        y0=y-0.3,
        y1=y,
        line_width=0,
        fillcolor=fillcolor,
        opacity=opacity,
    )

def _plot_gene_exon(fig, gene, gene_df, lfc, exon_color):

    """"""
    
    exon_plot_y = _get_gene_y_position(lfc)
    
    _start, _end = gene_df.Start.min(), gene_df.End.max()
    exon_plot_x = np.mean([float(_start), float(_end)])
    
    fig.add_annotation(x=exon_plot_x, y=exon_plot_y + 0.15, text=gene, showarrow=False,)
    fig.update_layout(showlegend=False)
    
    fig.add_shape(
        x0=_start,
        x1=_end,
        y0=exon_plot_y-0.18,
        y1=exon_plot_y-.12,
        line_width=0,
        fillcolor="red",
        opacity=1,
    )
    
    
    for exon in range(len(gene_df)):

        _plot_gene_single_exon(
            fig,
            gene_df.Start[exon],
            gene_df.End[exon],
            y=exon_plot_y,
            fillcolor=exon_color,
            opacity=0.25,
        )