"""
Author: Alexandra Lee
Date Created: 17 April 2020

Scripts to process expression data for pilot experiment
"""
import pandas as pd
import gzip
from Bio import SeqIO


def get_core_acc_genes(gene_annot_file):
    gene_annot = pd.read_csv(gene_annot_file, header=0, sep="\t", index_col=0)

    # Group genes by core and accessory annotation
    core_gene_ids = list(gene_annot[gene_annot["annotation"] == "core"].index)
    acc_gene_ids = list(gene_annot[gene_annot["annotation"] == "accessory"].index)

    return core_gene_ids, acc_gene_ids


def get_sample_grps(sample_annot_file):
    sample_annot = pd.read_csv(sample_annot_file, header=0, sep="\t", index_col=0)

    # Group genes by core and accessory annotation
    pao1_samples = list(sample_annot[sample_annot["genotype"] == "PAO1"].index)
    pa14_samples = list(sample_annot[sample_annot["genotype"] == "PA14"].index)

    return pao1_samples, pa14_samples


def dict_gene_num_to_ids(fasta_file):

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
