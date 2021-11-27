
import pandas as pd

def _gene_from_gtf(gtf_file, chromosome, gene):

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


def _get_feature(df, feature):
    return df.loc[df[2] == feature]


def _get_gene_body_bounds(gene_df):

    gene_body_min = gene_df[3].astype(int).min()
    gene_body_max = gene_df[4].astype(int).max()

    return gene_body_min, gene_body_max