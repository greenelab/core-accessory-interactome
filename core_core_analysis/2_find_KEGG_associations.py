# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1+dev
#   kernelspec:
#     display_name: Python [conda env:core_acc] *
#     language: python
#     name: conda-env-core_acc-py
# ---

# # Find KEGG associations
#
# This notebokk will create a table that has the KEGG pathways that are associated with the most stable and least stable core genes.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
from core_acc_modules import paths, utils, modules

random.seed(1)
# -

# Output files
pao1_core_stable_similarity_filename = "pao1_core_stable_associations.tsv"
pa14_core_stable_similarity_filename = "pa14_core_stable_associations.tsv"

# +
# Load transcriptional similarity df
pao1_similarity_scores_filename = "pao1_similarity_scores.tsv"
pa14_similarity_scores_filename = "pa14_similarity_scores.tsv"

pao1_similarity_scores = pd.read_csv(
    pao1_similarity_scores_filename, sep="\t", header=0, index_col=0
)
pa14_similarity_scores = pd.read_csv(
    pa14_similarity_scores_filename, sep="\t", header=0, index_col=0
)
# -

pao1_similarity_scores.head()

# +
# Load KEGG pathway data
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", index_col=0, header=None)
# -

pao1_pathways[2] = pao1_pathways[2].str.split(";").apply(set)
pao1_pathways.index = pao1_pathways.index.str.split(" - ").str[0]
pao1_pathways.head()

# ## Pathway annotations to PA14
#
# The annotations we have are only for PAO1 genes, so we will map PAO1 core genes to PA14 core genes to add annotations to PA14. This is possible since we are focused on only core genes, which have homologs between PAO1 and PA14

pao1_annotation_filename = paths.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(pao1_annotation_filename, "pao1")

gene_mapping_pao1 = gene_mapping_pao1["PA14_ID"].to_frame()

# ## Get pathway associations for most and least stable genes

# Get most and least stable core genes
most_stable_genes = pao1_similarity_scores[
    pao1_similarity_scores["label"] == "most stable"
].index
least_stable_genes = pao1_similarity_scores[
    pao1_similarity_scores["label"] == "least stable"
].index


def get_associated_pathways(genes_):
    rows = []
    for gene_id in genes_:
        pathway_bool = [
            gene_id in pao1_pathways.loc[pathway, 2] for pathway in pao1_pathways.index
        ]
        found_pathways = list(pao1_pathways[pathway_bool].index)
        rows.append({"gene id": gene_id, "pathways present": found_pathways})
    return pd.DataFrame(rows).set_index("gene id")


most_stable_associations = get_associated_pathways(most_stable_genes)
most_stable_associations.head()

least_stable_associations = get_associated_pathways(least_stable_genes)
least_stable_associations.head()

# +
# Concatenate dataframes
pao1_all_associations = pd.concat([most_stable_associations, least_stable_associations])
pao1_all_associations.head()

# TO DO: Rename index col
# -

# Map PA14 gene ids
pa14_all_associations = pao1_all_associations.merge(
    gene_mapping_pao1, left_index=True, right_index=True
)
pa14_all_associations.set_index("PA14_ID", inplace=True)
pa14_all_associations.head()

# Merge KEGG associations with transcriptional similarity information
pao1_all_associations = pao1_all_associations.merge(
    pao1_similarity_scores, left_index=True, right_index=True, how="left"
)
pa14_all_associations = pa14_all_associations.merge(
    pa14_similarity_scores, left_index=True, right_index=True, how="left"
)

pao1_all_associations.head()

# Save
pao1_all_associations.to_csv(pao1_core_stable_similarity_filename, sep="\t")
pa14_all_associations.to_csv(pa14_core_stable_similarity_filename, sep="\t")

# **Takeaway:**
#
# Based on the pathways associated with the most and least stable core genes, we find that the most stable core genes tend to be associated with pathways related to cellular maintenance including: protein transport systems, ribosomes, metabolism, type III, IV secretion system which mediates virulence.
#
# There are far fewer KEGG pathways that least stable core genes are found to be associated with. The least stable core genes are mostly associated with different types of metabolism.
#
# A google doc containing the most and least stable core genes and some information about them is [here](https://docs.google.com/spreadsheets/d/1SqEyBvutfbsOTo4afg9GiEzP32ZKplkN1a6MpAQBvZI/edit?usp=sharing).
