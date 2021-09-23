"""
Author: Alexandra Lee
Date Created: 17 April 2020

Scripts to process expression data for pilot experiment
"""
import pandas as pd
import random
import gzip
from Bio import SeqIO


def get_pao1_pa14_gene_map(gene_annotation_filename, reference_genotype):
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
    # Check that genotype is set correctly
    assert reference_genotype.lower() in [
        "pao1",
        "pa14",
    ], "Reference genotype string needs to be either pao1 or pa14"

    # Read gene annotation
    gene_mapping = pd.read_csv(gene_annotation_filename, sep=",", header=0)

    reference_id = "PAO1_ID" if reference_genotype.lower() == "pao1" else "PA14_ID"
    non_reference_id = "PA14_ID" if reference_genotype.lower() == "pao1" else "PAO1_ID"

    # Accessory are genes that don't have a homologous gene
    unmapped_genes = gene_mapping[gene_mapping[non_reference_id].isna()]
    acc_genes = list(unmapped_genes[reference_id])

    # Add label for core genes
    gene_mapping.loc[~gene_mapping[reference_id].isin(acc_genes), "annotation"] = "core"

    # Add column with number of genes mapped in case user
    # would like to only consider genes with 1-1 mapping
    gene_mapping["num_mapped_genes"] = (
        gene_mapping[non_reference_id].str.split(", ").str.len()
    )

    gene_mapping.set_index(reference_id, inplace=True)

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


def get_my_core_acc_genes(pao1_annotation_filename, pa14_annotation_filename, my_expression_pao1_df, my_expression_pa14_df):
    """
    Returns lists of core and accessory genes in your dataset
    using both PAO1 and PA14 reference genotypes

    pao1_annotation_file: str
        File containing mapping between PAO1 and PA14 gene ids downloaded from
        BACTOME annotation tool: https://pseudomonas-annotator.shinyapps.io/pa_annotator/

        Columns: PAO1_ID,Name,Product.Name,PA14_ID

    pa14_annotation_file: str
        File containing mapping between PA14 and PAO1 gene ids downloaded from
        BACTOME annotation tool: https://pseudomonas-annotator.shinyapps.io/pa_annotator/

        Columns: PA14_ID,Name,Product.Name,PAO1_ID

    my_expression_pao1_df: dataframe
        Gene expression dataframe aligned using PAO1 reference

    my_expression_pa14_df: dataframe
        Gene expression dataframe aligned using PA14 reference
    """
    # Get mapping between PAO1 and PA14 genes using PAO1 and PA14 references
    gene_mapping_pao1 = get_pao1_pa14_gene_map(pao1_annotation_filename, "pao1")
    gene_mapping_pa14 = get_pao1_pa14_gene_map(pa14_annotation_filename, "pa14")

    # Get core genes: genes that have a homolog between PAO1 and PA14
    core_pao1_genes, core_pa14_genes = get_core_genes(
        gene_mapping_pao1,
        gene_mapping_pa14,
        False
    )

    print(f"Number of PAO1 core genes: {len(core_pao1_genes)}")
    print(f"Number of PA14 core genes: {len(core_pa14_genes)}")

    # Select only core genes that are included in my dataset
    pao1_ref_genes = my_expression_pao1_df.columns
    my_core_pao1_genes = list(set(core_pao1_genes).intersection(pao1_ref_genes))

    print(f"Number of PAO1 core genes in my dataset: {len(my_core_pao1_genes)}")

    # Select only core genes that are included in my dataset
    pa14_ref_genes = my_expression_pa14_df.columns
    my_core_pa14_genes = list(set(core_pa14_genes).intersection(pa14_ref_genes))

    print(f"Number of PA14 core genes in my dataset: {len(my_core_pa14_genes)}")

    # Get PAO1-specific genes
    pao1_acc = list(set(pao1_ref_genes) - set(my_core_pao1_genes))
    print(f"Number of PAO1-specific genes: {len(pao1_acc)}")

    # Get PA14-specific genes
    pa14_acc = list(set(pa14_ref_genes) - set(my_core_pa14_genes))
    print(f"Number of PA14-specific genes: {len(pa14_acc)}")

    return {"core_pao1": my_core_pao1_genes,
            "core_pa14": my_core_pa14_genes,
            "acc_pao1": pao1_acc,
            "acc_pa14": pa14_acc
            }


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
    for downstream analysis.

    For each sample in `input_data`, shuffle gene values (i.e. shuffle values independently per sample)
    in order to disrupt any gene structure (i.e. underlying relationship between genes).

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


def check_sample_ordering(expression_file, metadata_file):
    """
    This function checks that the ordering of the samples matches
    between the expression file and the metadata file. This
    ordering is used for calculating DEGs.
    """
    # Check ordering of sample ids is consistent between gene expression data and metadata
    metadata = pd.read_csv(metadata_file, sep="\t", header=0, index_col=0)
    metadata_sample_ids = list(metadata.index)

    expression_data = pd.read_csv(expression_file, sep="\t", header=0, index_col=0)
    expression_sample_ids = list(expression_data.index)

    if metadata_sample_ids == expression_sample_ids:
        print("sample ids are ordered correctly")
    else:
        # Convert gene expression ordering to be the same as
        # metadata sample ordering
        print("sample ids don't match, going to re-order gene expression samples")
        expression_data = expression_data.loc[metadata_sample_ids]
        expression_data.to_csv(expression_file, sep="\t")


def get_sample_ids(metadata_filename, experiment_colname, sample_colname, experiment_id):
    """
    Returns sample ids (found in gene expression df) associated with
    a given list of experiment ids (found in the metadata)

    Arguments
    ----------
    metadata_filename: str
        File containing metadata
    experiment_colname: str
        Column header that contains experiment id that maps expression data
        and metadata
    sample_colname: str
        Column header that contains sample id that maps expression data
        and metadata
    experiment_id: str
        Selected experiment id to grab samples from

    """
    # Read in metadata
    metadata = pd.read_csv(metadata_filename, header=0)
    metadata.set_index(experiment_colname, inplace=True)

    selected_metadata = metadata.loc[experiment_id]
    sample_ids = list(selected_metadata[sample_colname])

    return sample_ids
