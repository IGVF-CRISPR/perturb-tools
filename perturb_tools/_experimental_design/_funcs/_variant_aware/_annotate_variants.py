
import pandas as pd
import glob
import vintools as v

def _annotate_variants_single_cell_line(df, path_to_vcf, cell_line_name):

    """
    Adds three columns to input df:
        1. column indicating boolean presence of a variant in the sgRNA sequence
        2. column indicating reference allele
        3. column indicating variant allele detected
    """

    vcf = v.ut.read_vcf(path_to_vcf)

    variant_pos_all_chromosomes = []
    all_chromosome_sgRNA_bins = []

    # loop through chromosomes within the list of sgRNAs
    for chromosome in df.chr.unique():
        chromosome_df = df.loc[df.chr == chromosome]

        # get sgRNA "bins" on a single chromosome
        individual_bins = []

        # loop through guides wtihin a chromosome
        for guide in range(len(chromosome_df)):
            left = chromosome_df.guide_begin.iloc[guide]
            right = chromosome_df.guide_end.iloc[guide]

            individual_bins.append(np.array([left, right]))

        # filter variant locations into sgRNA "bins"
        single_chromosome_sgRNA_bins = var_in_guide_ranges(vcf.POS, individual_bins)
        all_chromosome_sgRNA_bins.append(single_chromosome_sgRNA_bins)

        variant_pos = []
        variant_ref = []
        variant_alt = []

        for guide_bin in single_chromosome_sgRNA_bins:
            if guide_bin.sum() == 0:
                variant_pos.append(False)
                variant_ref.append(False)
                variant_alt.append(False)
            else:
                # annotate boolean presence
                variant_pos.append(True)
                # annotate reference allele
                variant_ref.append(vcf.iloc[np.where(guide_bin == True)]["REF"].values)
                # annotate variant allele
                variant_alt.append(vcf.iloc[np.where(guide_bin == True)]["ALT"].values)
        #                 variant_tally_single_chromosome.append(
        #                     vcf.iloc[np.where(guide_bin == True)]["Cell_line"].values[0]
        #                 )

        # add a column titled by the name of the cell line; contains boolean variant information

        variant_pos_col = "_".join([cell_line_name, "_variant_POS"])
        variant_ref_col = "_".join([cell_line_name, "_variant_REF"])
        variant_alt_col = "_".join([cell_line_name, "_variant_ALT"])

        chromosome_df[variant_pos_col] = variant_pos
        chromosome_df[variant_ref_col] = variant_ref
        chromosome_df[variant_alt_col] = variant_alt

        variant_pos_all_chromosomes.append(chromosome_df)

    df = pd.concat(variant_pos_all_chromosomes)

    return df


def _annotate_variants(df, vcf_dir):
    
    """Wraps single cell line function for each cell line."""
    
    vcf_dir_full_path = "".join([vcf_dir, "*.vcf"])
    vcf_path_list = glob.glob(vcf_dir_full_path)

    all_cell_lines_dfs = []

    for vcf_file in vcf_path_list:
        name = vcf_file.split(".vcf")[-2].split("/")[-1]
        print("Now annotating variants from", name, "...")
        single_cell_line_df = _annotate_variants_single_cell_line(
            df, vcf_file, cell_line_name=name
        )
        all_cell_lines_dfs.append(single_cell_line_df)

    df = pd.concat(all_cell_lines_dfs).drop_duplicates("protospacer").fillna(False)

    return df