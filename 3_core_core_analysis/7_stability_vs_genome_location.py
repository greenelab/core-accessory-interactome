# -*- coding: utf-8 -*-
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

# # Stability vs genome location
#
# Here we want to know about the relationship between stability and genome distance: Are more stable genes located in similar locations on the genome in PAO1 vs PA14 compared to least stable genes?
#
# We hypothesize that least stable genes are not syntenic (i.e. located in different regions of the genome) across strains, which might indicate a different some transcriptional re-wiring.
#
# There are 2 approaches we take here:
# 1. For each least/most stable gene, get the neighboring core/homologous genes and determine if they match between PAO1 vs PA14.
# 2. Scale PA14 gene ids to same range as PAO1 then take distance between PAO1 gene and homolog of PAO1 gene
#
# https://academic.oup.com/bioinformatics/article/36/Supplement_1/i21/5870523
# https://pubmed.ncbi.nlm.nih.gov/23323735/
# https://arxiv.org/pdf/1307.4291.pdf

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import paths, utils

random.seed(1)
# -

# User param
window = 2

# ## Load data

# Input similarity scores and annotations filenames
# Since the results are similar we only need to look at the scores for one strain type
pao1_similarity_filename = "pao1_core_similarity_associations_final_spell.tsv"
pa14_similarity_filename = "pa14_core_similarity_associations_final_spell.tsv"

# Import df
pao1_similarity = pd.read_csv(pao1_similarity_filename, sep="\t", index_col=0, header=0)
pa14_similarity = pd.read_csv(pa14_similarity_filename, sep="\t", index_col=0, header=0)

print(pao1_similarity.shape)
pao1_similarity.head()

print(pa14_similarity.shape)
pa14_similarity.head()

# ## Get most and least stable genes

pao1_least_stable = pao1_similarity[pao1_similarity["label"] == "least stable"].index
pao1_most_stable = pao1_similarity[pao1_similarity["label"] == "most stable"].index

# ## Get mapping from PAO1 to PA14 ids

pao1_annotation_filename = paths.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(pao1_annotation_filename, "pao1")

gene_mapping_pao1 = gene_mapping_pao1["PA14_ID"].to_frame()

gene_mapping_pao1.head()


# ## Approach 1
#
# For least stable gene:
# * select window of core genes
# * compare sets

# TO DO:
# Which all core set to use
# Do this for a random set of core genes as a baseline
def test(
    pao1_core_genes_df,
    pa14_core_genes_df,
    gene_id_mapping_df,
    window_size,
    core_genes_to_examine,
):
    """
    Function that compares ...

    gene_id_mapping_df: df
        Mapping PAO1 ids to PA14 ids

    core_genes_to_examine: df
        Least or most stable core genes using PAO1 ids
    """
    # Sort to make sure core ids are in order
    pao1_core_genes_df = pao1_core_genes_df.sort_index()
    pa14_core_genes_df = pa14_core_genes_df.sort_index()

    # List of percent overlap
    percent_overlap_lst = []
    for pao1_id in core_genes_to_examine:
        print(pao1_id)

        # Get PAO1 neighboring core genes
        pao1_gene_idx = pao1_core_genes_df.index.get_loc(pao1_id)
        pao1_gene_neighborhood_ids = pao1_core_genes_df.iloc[
            pao1_gene_idx - window_size : pao1_gene_idx + window_size + 1
        ].index

        # Remove least/most stable gene
        pao1_gene_neighborhood_ids = pao1_gene_neighborhood_ids.drop(labels=pao1_id)

        print(pao1_gene_idx)
        print(pao1_gene_neighborhood_ids)

        # Map PAO1 least/most stable gene to PA14 id
        mapped_pa14_id = gene_id_mapping_df.loc[pao1_id, "PA14_ID"]
        print(mapped_pa14_id)

        # Get PA14 neighboring core genes
        pa14_gene_idx = pa14_core_genes_df.index.get_loc(mapped_pa14_id)
        pa14_gene_neighborhood_ids = pa14_core_genes_df.iloc[
            pa14_gene_idx - window_size : pa14_gene_idx + window_size + 1
        ].index

        # Remove least/most stable gene
        pa14_gene_neighborhood_ids = pa14_gene_neighborhood_ids.drop(
            labels=mapped_pa14_id
        )

        print(pa14_gene_idx)
        print(pa14_gene_neighborhood_ids)

        # Convert neighboring PAO1 ids to PA14 ids
        mapped_pao1_gene_neighborhood_ids = gene_id_mapping_df.loc[
            pao1_gene_neighborhood_ids, "PA14_ID"
        ]
        print(mapped_pao1_gene_neighborhood_ids)

        # Compare gene ids
        overlap_neighborhood_ids = set(pa14_gene_neighborhood_ids).intersection(
            set(mapped_pao1_gene_neighborhood_ids)
        )
        print("overlap", overlap_neighborhood_ids)

        # Save percent matched core genes
        percent_overlap = len(overlap_neighborhood_ids) / (2 * window_size)
        print(len(overlap_neighborhood_ids))
        percent_overlap_lst.append(percent_overlap)

    return percent_overlap_lst


least_matched_neighborhood = test(
    pao1_similarity, pa14_similarity, gene_mapping_pao1, window, pao1_least_stable
)

most_matched_neighborhood = test(
    pao1_similarity, pa14_similarity, gene_mapping_pao1, window, pao1_most_stable
)

# Random set of genes
random_gene_set = pao1_similarity.sample(len(pao1_least_stable)).index
random_matched_neighborhood = test(
    pao1_similarity, pa14_similarity, gene_mapping_pao1, window, random_gene_set
)

least_matched_neighborhood

most_matched_neighborhood


# ## Approach 2

# +
# Approach 2
# Calculate distance for all least/most genes
# Compare sets
# -

def test2(
    pao1_core_genes_df,
    pa14_core_genes_df,
    gene_id_mapping_df,
    window_size,
    core_genes_to_examine,
):
    """
    Function that compares ...

    gene_id_mapping_df: df
        Mapping PAO1 ids to PA14 ids

    core_genes_to_examine: df
        Least or most stable core genes using PAO1 ids
    """
    diff_lst = []

    # Sort to make sure core ids are in order
    pao1_core_genes_df = pao1_core_genes_df.sort_index()
    pa14_core_genes_df = pa14_core_genes_df.sort_index()

    # Get number corresponding to PA id and scale to be between 0 - 1
    pao1_core_genes_df["PAO1_num"] = (
        pao1_core_genes_df.index.str.split("PA").str[-1].astype("float")
    )
    pao1_core_genes_df["PAO1_scaled"] = pao1_core_genes_df["PAO1_num"] / max(
        pao1_core_genes_df["PAO1_num"]
    )

    pa14_core_genes_df["PA14_num"] = (
        pa14_core_genes_df.index.str.split("PA14_").str[-1].astype("float")
    )
    pa14_core_genes_df["PA14_scaled"] = pa14_core_genes_df["PA14_num"] / max(
        pa14_core_genes_df["PA14_num"]
    )

    print(pao1_core_genes_df.head())
    print(pao1_core_genes_df.tail())

    print(pa14_core_genes_df.head())
    print(pa14_core_genes_df.tail())

    # Calculate the distance between most/least stable genes in PAO1 vs PA14
    # using scaled location
    for pao1_id in core_genes_to_examine:
        print(pao1_id)

        pao1_loc = pao1_core_genes_df.loc[pao1_id, "PAO1_scaled"]

        # Map PAO1 least/most stable gene to PA14 id
        mapped_pa14_id = gene_id_mapping_df.loc[pao1_id, "PA14_ID"]
        print(mapped_pa14_id)

        pa14_loc = pa14_core_genes_df.loc[mapped_pa14_id, "PA14_scaled"]

        # Calculate the difference
        pao1_vs_pa14_loc = abs(pao1_loc - pa14_loc)

        diff_lst.append(pao1_vs_pa14_loc)

    return diff_lst


pao1_least_dist = test2(
    pao1_similarity, pa14_similarity, gene_mapping_pao1, window, pao1_least_stable
)

pao1_most_dist = test2(
    pao1_similarity, pa14_similarity, gene_mapping_pao1, window, pao1_most_stable
)

pao1_random_dist = test2(
    pao1_similarity, pa14_similarity, gene_mapping_pao1, window, random_gene_set
)

# ## Plot

# +
plt.figure(figsize=(10, 8))
sns.violinplot(data=least_matched_neighborhood, color="grey")
sns.violinplot(data=most_matched_neighborhood, color="blue")
sns.violinplot(data=random_matched_neighborhood, color="pink")

plt.title("Stability of secretion system genes", fontsize=14)
plt.legend()

# +
plt.figure(figsize=(10, 8))
sns.violinplot(data=pao1_least_dist, color="grey")
sns.violinplot(data=pao1_most_dist, color="blue")
sns.violinplot(data=pao1_random_dist, color="pink")

plt.title("Stability of secretion system genes", fontsize=14)
plt.legend()

# +
# TO DO
# Check which core genes to start with
# Check calculations for Approach 1, control?
# Check calculations for Approach 2, control?
# Fix plotting
# Synteny calculations lookup
