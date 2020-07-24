"""
Author: Alexandra Lee
Date Created: 17 April 2020

Scripts to process expression data for pilot experiment
"""
import pandas as pd
import random
import gzip
from Bio import SeqIO

random.seed(123)


def get_pao1_pa14_gene_map(gene_annotation_file, reference_genotype):
    """
    Returns file mapping PAO1 gene ids to PA14 ids and label which genes are core and
    the number of genes that are mapped

    Arguments
    ----------
    gene_annotation_file: str
        File containing mapping between PAO1 and PA14 gene ids downloaded from
        BACTOME annotation tool: https://pseudomonas-annotator.shinyapps.io/pa_annotator/

        Columns: PAO1_ID,Name,Product.Name,PA14_ID

    reference_genotype: str
        Either 'pao1' or 'pa14'

    """

    # Read gene annotation
    gene_mapping = pd.read_csv(gene_annotation_file, sep=",", header=0)

    # Accessory are genes that don't have a mapping to the PA14 ID
    if reference_genotype.lower() == "pao1":
        unmapped_genes = gene_mapping[gene_mapping["PA14_ID"].isna()]
        acc_genes = list(unmapped_genes["PAO1_ID"])

        # Add label for core genes
        gene_mapping.loc[
            ~gene_mapping["PAO1_ID"].isin(acc_genes), "annotation"
        ] = "core"

        # Add column with number of genes mapped in case user
        # would like to only consider genes with 1-1 mapping
        gene_mapping["num_mapped_genes"] = (
            gene_mapping["PA14_ID"].str.split(", ").str.len()
        )

        gene_mapping.set_index("PAO1_ID", inplace=True)

    elif reference_genotype.lower() == "pa14":
        unmapped_genes = gene_mapping[gene_mapping["PAO1_ID"].isna()]
        acc_genes = list(unmapped_genes["PA14_ID"])

        # Add label for core genes
        gene_mapping.loc[
            ~gene_mapping["PA14_ID"].isin(acc_genes), "annotation"
        ] = "core"

        # Add column with number of genes mapped in case user
        # would like to only consider genes with 1-1 mapping
        gene_mapping["num_mapped_genes"] = (
            gene_mapping["PAO1_ID"].str.split(", ").str.len()
        )

        gene_mapping.set_index("PA14_ID", inplace=True)

    return gene_mapping


def get_core_genes(pao1_ref_mapping_df, pa14_ref_mapping_df, is_mapping_1_1):
    """
    Returns list of core genes using PAO1 ids and PA14 ids

    Arguments
    ----------
    pao1_ref_mapping_df: df
        Dataframe generated from get_pao1_pa14_gene_map() to give mapping
        from PAO1 ids to PA14 ids

        Columns: PAO1_ID, Name, Product.Name, PA14_ID, annotation, num_mapped_genes

    pa14_ref_mapping_df: df
        Dataframe generated from get_pao1_pa14_gene_map() to give mapping
        from PA14 ids to PAO1 ids

        Columns: PA14_ID, Name, Product.Name, PAO1_ID, annotation, num_mapped_genes

    is_mapping_1_1: bool
        True if only want to return core genes that have a 1-1 mapping between PAO1 and PA14
    """

    if is_mapping_1_1:
        # Include only core genes that have a 1-1 gene mapping
        pao1_ref_core_df = pao1_ref_mapping_df[
            (pao1_ref_mapping_df["annotation"] == "core")
            & (pao1_ref_mapping_df["num_mapped_genes"] == 1.0)
        ]

        pa14_ref_core_df = pa14_ref_mapping_df[
            (pa14_ref_mapping_df["annotation"] == "core")
            & (pa14_ref_mapping_df["num_mapped_genes"] == 1.0)
        ]
    else:
        pao1_ref_core_df = pao1_ref_mapping_df[
            (pao1_ref_mapping_df["annotation"] == "core")
        ]
        pa14_ref_core_df = pa14_ref_mapping_df[
            (pa14_ref_mapping_df["annotation"] == "core")
        ]

    # Get list of pao1 core genes
    pao1_core_genes = pd.DataFrame(
        list(pao1_ref_core_df.index) + list(pa14_ref_core_df["PAO1_ID"].values),
        columns=["gene id"],
    )
    # Reshape to get single value per row
    pao1_core_genes = pd.DataFrame(
        pao1_core_genes["gene id"].str.split(", ").sum(), columns=["gene id"]
    )
    # Remove duplicates that might exist after taking the union
    pao1_core_genes.drop_duplicates(keep="first", inplace=True)
    pao1_core_genes = list(pao1_core_genes["gene id"])

    # Get list of pa14 core genes
    pa14_core_genes = pd.DataFrame(
        list(pa14_ref_core_df.index) + list(pao1_ref_core_df["PA14_ID"].values),
        columns=["gene id"],
    )
    # Reshape to get single value per row
    pa14_core_genes = pd.DataFrame(
        pa14_core_genes["gene id"].str.split(", ").sum(), columns=["gene id"]
    )
    # Remove duplicates that might exist after taking the union
    pa14_core_genes.drop_duplicates(keep="first", inplace=True)
    pa14_core_genes = list(pa14_core_genes["gene id"])

    return pao1_core_genes, pa14_core_genes


def get_sample_grps(sample_annot_file):
    """
    Returns list of sample ids that are PAO1 and PA14 samples based on
    the experiment metadata

    Arguments
    ----------
    sample_annot_file: str
        File containing two columns: accession, genotype
    """

    sample_annot = pd.read_csv(sample_annot_file, header=0, sep="\t", index_col=0)

    # Group genes by core and accessory annotation
    pao1_samples = list(sample_annot[sample_annot["genotype"] == "PAO1"].index)
    pa14_samples = list(sample_annot[sample_annot["genotype"] == "PA14"].index)

    return pao1_samples, pa14_samples


def dict_gene_num_to_ids(fasta_file):
    """
    Returns mapping between gene number created from SRA (i.e. PGD######)
    and gene locus (PA######). The PA###### is the position of the gene on
    the Pseudomonas genome and is more interpretable compared to the
    SRA-generated PGD######

    Arguments
    ----------
    fasta_file: str
        Reference transcriptome file 
    """

    seq_id_to_gene_id = {}

    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id
            seq_record_desc = record.description
            seq_record_desc_lst = seq_record_desc.split(";")
            tagged_item = [
                item.split("=")[1]
                for item in seq_record_desc_lst
                if "locus_tag" in item
            ][0]

            seq_id_to_gene_id[seq_id] = tagged_item

    return seq_id_to_gene_id


def permute_expression_data(input_data):
    """
    Returns permuted version of input expression data to be used as a baseline
    for downstream analysis

    Arguments
    ----------
    input_data: df
        File containing normalized gene expression data

        ------------------------------| PA0001 | PA0002 |...
        05_PA14000-4-2_5-10-07_S2.CEL | 0.8533 | 0.7252 |...
        54375-4-05.CEL                | 0.7789 | 0.7678 |...
        ...                           | ...    | ...    |...

    Returns
    -------
        Permuted expression data
    """

    # For each sample, shuffle gene values (i.e. shuffle values independently per sample) in order to
    # disrupt any gene structure
    shuffled_arr = []
    num_samples = input_data.shape[0]

    for i in range(num_samples):
        row = list(input_data.values[i])
        shuffled_row = random.sample(row, len(row))
        shuffled_arr.append(shuffled_row)

    shuffled_data = pd.DataFrame(
        shuffled_arr, index=input_data.index, columns=input_data.columns
    )

    return shuffled_data
