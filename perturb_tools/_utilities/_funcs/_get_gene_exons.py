
import pandas as pd

from ..._experimental_design._funcs._GTF_Module import _GTF

def _get_gene_df_from_gtf(gtf_file, gene, chromosome):

    """"""

    GeneDict = {}

    n = 0
    for i, line in enumerate(gtf_file.readlines()):
        line = line.strip("\n").split("\t")
        if line[0] == chromosome:
            line_ = line[8].split(" ")
            if line_[10] == "gene_name":
                if line_[11].startswith('"{}'.format(gene)):
                    GeneDict[n] = line
                    n += 1
    df = pd.DataFrame.from_dict(GeneDict, orient="index")

    return df


def _get_feature(df, feature="exon"):

    """"""

    df_ = (
        df.loc[df[2] == feature]
        .loc[df[1] == "HAVANA"]
        .reset_index(drop=True)[[0, 3, 4]]
    )
    df_.columns = ["Chromosome", "Start", "End"]
    return df_


def _get_gene_exons(gene, chromosome, gtf_path):

    """"""

    gtf = _GTF(gtf_path)
    gene_df = _get_feature(_get_gene_df_from_gtf(gtf.file, gene, chromosome))

    return gene_df