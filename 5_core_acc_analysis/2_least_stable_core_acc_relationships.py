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

# # Accessory genes related to least stable core genes
#
# Since we found that least stable core genes are more co-expressed with accessory genes. Let's look at who those accessory genes are. In our [previous notebook](../3_core_core_analysis/4_find_related_acc_genes.ipynb) we annotated least and most stable core genes with their top co-expressed accessory gene.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import scipy
import pandas as pd
import numpy as np
from scripts import paths, utils

random.seed(1)
# -

# Load transcriptional similarity df
# These are the subset of genes that we will consider
pao1_similarity_scores_filename = (
    "../3_core_core_analysis/pao1_core_similarity_associations_final_spell.tsv"
)
pa14_similarity_scores_filename = (
    "../3_core_core_analysis/pa14_core_similarity_associations_final_spell.tsv"
)

pao1_similarity_scores = pd.read_csv(
    pao1_similarity_scores_filename, sep="\t", header=0, index_col=0
)
pa14_similarity_scores = pd.read_csv(
    pa14_similarity_scores_filename, sep="\t", header=0, index_col=0
)

# Get least stable core genes
pao1_least_stable_genes = list(
    pao1_similarity_scores[pao1_similarity_scores["label"] == "least stable"].index
)
pa14_least_stable_genes = list(
    pa14_similarity_scores[pa14_similarity_scores["label"] == "least stable"].index
)

# The co-expressed accessory genes are listed in the `Related acc genes` column. The values are a list of accessory genes that were in the top 10 co-expressed genes, otherwise "No accessory genes"

pao1_least_df = pao1_similarity_scores.loc[pao1_least_stable_genes]

pa14_least_df = pa14_similarity_scores.loc[pa14_least_stable_genes]

pao1_least_df.head()

# ## Concatenate accessory genes into a list
#
# Note: We have to use `eval` here because we have a mix of strings and lists in our column of interest. In the future we could use an empty list instead of a string.

pao1_least_processed_df = pao1_least_df[
    pao1_least_df["Related acc genes"] != "No accessory genes"
]["Related acc genes"]
pa14_least_processed_df = pa14_least_df[
    pa14_least_df["Related acc genes"] != "No accessory genes"
]["Related acc genes"]

pao1_least_counts_df = (
    pd.Series(pao1_least_processed_df.apply(eval).sum())
    .value_counts()
    .to_frame("counts")
)
pa14_least_counts_df = (
    pd.Series(pa14_least_processed_df.apply(eval).sum())
    .value_counts()
    .to_frame("counts")
)

pao1_least_counts_df

pa14_least_counts_df

# ## Add gene names

pao1_annotation_filename = paths.GENE_PAO1_ANNOT
pa14_annotation_filename = paths.GENE_PA14_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(pao1_annotation_filename, "pao1")
gene_mapping_pa14 = utils.get_pao1_pa14_gene_map(pa14_annotation_filename, "pa14")

pao1_gene_name_map = gene_mapping_pao1["Name"].to_frame()
pa14_gene_name_map = gene_mapping_pa14["Name"].to_frame()

gene_mapping_pao1.head()

# Add gene names
pao1_least_counts_df = pao1_least_counts_df.merge(
    gene_mapping_pao1["Name"], left_index=True, right_index=True, how="left"
)
pa14_least_counts_df = pa14_least_counts_df.merge(
    pa14_gene_name_map["Name"], left_index=True, right_index=True, how="left"
)

pao1_least_counts_df

pa14_least_counts_df

# Save
pao1_least_counts_df.to_csv("pao1_acc_coexpressed_with_least_stable.tsv", sep="\t")
pa14_least_counts_df.to_csv("pa14_acc_coexpressed_with_least_stable.tsv", sep="\t")
