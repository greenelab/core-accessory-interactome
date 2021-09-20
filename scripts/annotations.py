"""
Author: Alexandra Lee
Date Created: 20 September 2021

Scripts to annotate genes
"""
import pandas as pd
import utils


def load_format_operons(operon_filename):
    """
    This function inputs operon annotations from the `operon_filename`
    and returns a pandas df with gene ids as index and the associated
    operon as the column values.
    """

    operon_df = pd.read_csv(operon_filename, index_col=0, header=0)

    # There are 247 PAO1 genes with multiple annotations
    # This operon df contains annotations from predicted operons based on DOOR database
    # predictions which make up the majority of the operons) as well as some that
    # are curated (i.e. PseudoCAP)
    # There are some that have multiple PseudoCAP annotations too

    # Here we will keep the last PseudoCAP annotations
    # To ensure that the PseudoCAP annotations are the last ones, we will sort the values
    operon_df = operon_df.sort_values(by=["locus_tag", "source_database"])
    operon_df = operon_df.set_index("locus_tag")
    operon_df = operon_df[~operon_df.index.duplicated(keep="last")]

    # Only include columns for gene id and operon_name
    operon_final_df = operon_df["operon_name"].to_frame()

    return operon_final_df


def load_format_KEGG(kegg_filename):
    """
    This function inputs operon annotations from the `kegg_filename`
    and returns a pandas df with KEGG pathways as index and the associated
    genes as the column values.
    """

    kegg_df = pd.read_csv(kegg_filename, index_col=0, header=None)

    kegg_df[2] = kegg_df[2].str.split(";").apply(set)
    kegg_df.index = kegg_df.index.str.split(" - ").str[0]

    return kegg_df


def map_core_acc_annot(
    pao1_membership_df,
    pa14_membership_df,
    pao1_expression_filename,
    pa14_expression_filename,
    pao1_annot_filename,
    pa14_annot_filename,
):
    """
    This returns a pandas df with gene ids as index and label
    for whether that gene is "core" or "acc"
    """

    pao1_expression = pd.read_csv(
        pao1_expression_filename, sep="\t", index_col=0, header=0
    )
    pa14_expression = pd.read_csv(
        pa14_expression_filename, sep="\t", index_col=0, header=0
    )
    core_acc_dict = utils.get_my_core_acc_genes(
        pao1_annot_filename, pa14_annot_filename, pao1_expression, pa14_expression
    )

    pao1_core = core_acc_dict["core_pao1"]
    pa14_core = core_acc_dict["core_pa14"]
    pao1_acc = core_acc_dict["acc_pao1"]
    pa14_acc = core_acc_dict["acc_pa14"]

    pao1_membership_df.loc[pao1_core, "core/acc"] = "core"
    pao1_membership_df.loc[pao1_acc, "core/acc"] = "acc"
    pa14_membership_df.loc[pa14_core, "core/acc"] = "core"
    pa14_membership_df.loc[pa14_acc, "core/acc"] = "acc"

    # Make sure to sort by gene id
    # NOTE PA14 gene ids don't increment by 1, but by 10 or 20 are we missing some genes?
    pao1_arr = pao1_membership_df.sort_index()
    pa14_arr = pa14_membership_df.sort_index()

    return pao1_arr, pa14_arr
