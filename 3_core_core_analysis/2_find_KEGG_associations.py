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
# This notebook will create a table that has the KEGG pathways that are associated with all genes, but we are particularly interested in those that are associated with the most and least stable genes.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
from scripts import paths, utils, modules, annotations

random.seed(1)
# -

# Output files
pao1_out_filename = "pao1_core_similarity_associations.tsv"
pa14_out_filename = "pa14_core_similarity_associations.tsv"

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

print(pao1_similarity_scores.shape)
pao1_similarity_scores.head()

print(pa14_similarity_scores.shape)
pa14_similarity_scores.head()

# Load KEGG pathway data
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

pao1_pathways = annotations.load_format_KEGG(pao1_pathway_filename)
pao1_pathways.head()

# ## Pathway annotations to PA14
#
# The annotations we have are only for PAO1 genes, so we will map PAO1 core genes to PA14 core genes to add annotations to PA14. This is possible since we are focused on only core genes, which have homologs between PAO1 and PA14

pao1_annotation_filename = paths.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(pao1_annotation_filename, "pao1")

gene_mapping_pao1 = gene_mapping_pao1["PA14_ID"].to_frame()


# ## Get pathway associations for all genes

def get_associated_pathways(genes_):
    rows = []
    for gene_id in genes_:
        pathway_bool = [
            gene_id in pao1_pathways.loc[pathway, 2] for pathway in pao1_pathways.index
        ]
        found_pathways = list(pao1_pathways[pathway_bool].index)
        rows.append({"gene id": gene_id, "pathways present": found_pathways})
    return pd.DataFrame(rows).set_index("gene id")


# Get KEGG associations for all genes in PAO1
all_pao1_gene_ids = list(pao1_similarity_scores.index)
pao1_associations = get_associated_pathways(all_pao1_gene_ids)

print(pao1_associations.shape)
pao1_associations.head()

# Map PA14 gene ids
pa14_associations = pao1_associations.merge(
    gene_mapping_pao1, left_index=True, right_index=True
)
pa14_associations.set_index("PA14_ID", inplace=True)
print(pa14_associations.shape)
pa14_associations.head()

# Merge KEGG associations with transcriptional similarity information
pao1_associations = pao1_similarity_scores.merge(
    pao1_associations, left_index=True, right_index=True, how="left"
)
pa14_associations = pa14_similarity_scores.merge(
    pa14_associations, left_index=True, right_index=True, how="left"
)

print(pao1_associations.shape)
pao1_associations.head()

print(pa14_associations.shape)
pa14_associations.head()

# Save
pao1_associations.to_csv(pao1_out_filename, sep="\t")
pa14_associations.to_csv(pa14_out_filename, sep="\t")

# **Takeaway:**
#
# Based on the pathways associated with the most and least stable core genes, we find that the most stable core genes tend to be associated with pathways related to cellular maintenance including: protein transport systems, ribosomes, metabolism, type III, IV secretion system which mediates virulence.
#
# There are far fewer KEGG pathways that least stable core genes are found to be associated with. The least stable core genes are mostly associated with different types of metabolism.
#
# A google doc containing the most and least stable core genes and some information about them is [here](https://docs.google.com/spreadsheets/d/1SqEyBvutfbsOTo4afg9GiEzP32ZKplkN1a6MpAQBvZI/edit?usp=sharing).
