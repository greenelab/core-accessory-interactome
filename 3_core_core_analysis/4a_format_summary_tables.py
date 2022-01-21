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

# # Format summary tables
#
# This notebook inputs the summary tables for PAO1 and PA14 and formats it to be one table using transcriptional similarity values for PAO1-aligned genes and includes expression statistics for both PAO1-aligned and PA14-aligned genes

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
from scripts import paths, utils

random.seed(1)
# -

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

# Select expression statistics for PA14-aligned data
pa14_subset = pa14_similarity[
    [
        "PAO1 homolog id",
        "mean expression",
        "standard deviation expression",
        "25% expression",
        "50% expression",
        "75% expression",
        "min expression",
        "max expression",
        "variance expression",
        "range expression",
        "pathways present",
        "Related acc genes",
    ]
].set_index("PAO1 homolog id")

print(pa14_subset.shape)
pa14_subset.head()

# Merge dataframes based on PAO1 gene id
all_similarity = pao1_similarity.merge(
    pa14_subset, left_index=True, right_index=True, suffixes=["_pao1", "_p14"]
)

print(all_similarity.shape)
all_similarity.head()

# +
# Some manual checks
# all_similarity.loc[["PA0744", "PA2588"]]
# pao1_similarity.loc[["PA0744", "PA2588"]]
# -

# Output
all_similarity.to_csv("all_core_similarity_associations_final_spell.tsv", sep="\t")
