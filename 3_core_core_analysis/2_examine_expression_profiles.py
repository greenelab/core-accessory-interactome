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

# # Examine transcription profiles
#
# This notebook tries to examine why genes are found to be "most stable" and "least stable." This notebook also performs a small exploratory analysis to check that the genes we are calling "most stable" and "least stable" are _real_. If genes are "most stable" because they are always on then this is not as interesting. To examine this we will add statistics about the expression distribution to the transcriptional similarity scores matrix.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import numpy as np
import pandas as pd
import seaborn as sns
from scripts import paths, utils, modules, annotations

random.seed(1)
# -

# Output files
pao1_out_filename = "pao1_core_similarity_expression_stats.tsv"
pa14_out_filename = "pa14_core_similarity_expression_stats.tsv"

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

# ### Add expression statistics

# Load expression data
pao1_expression_filename = paths.PAO1_COMPENDIUM
pa14_expression_filename = paths.PA14_COMPENDIUM

# Transpose gene expression matrices to be gene x sample
pao1_expression = pd.read_csv(
    pao1_expression_filename, sep="\t", header=0, index_col=0
).T
pa14_expression = pd.read_csv(
    pa14_expression_filename, sep="\t", header=0, index_col=0
).T

pao1_expression.head()

# Calculate mean, median, variance, range of expression
pao1_expression["median expression"] = pao1_expression.median(axis=1)
pao1_expression["mean expression"] = pao1_expression.mean(axis=1)
pao1_expression["variance expression"] = pao1_expression.var(axis=1)
pao1_expression["max expression"] = pao1_expression.max(axis=1)
pao1_expression["min expression"] = pao1_expression.min(axis=1)
pao1_expression["range expression"] = pao1_expression.max(axis=1) - pao1_expression.min(
    axis=1
)

pa14_expression["median expression"] = pa14_expression.median(axis=1)
pa14_expression["mean expression"] = pa14_expression.mean(axis=1)
pa14_expression["variance expression"] = pa14_expression.var(axis=1)
pa14_expression["max expression"] = pa14_expression.max(axis=1)
pa14_expression["min expression"] = pa14_expression.min(axis=1)
pa14_expression["range expression"] = pa14_expression.max(axis=1) - pa14_expression.min(
    axis=1
)

select_columns = pao1_expression.columns[
    pao1_expression.columns.str.contains("expression")
]

pao1_expression_stats = pao1_expression[select_columns]
pa14_expression_stats = pa14_expression[select_columns]

pao1_expression_stats

# Merge expression statistics with transcriptional similarity information
pao1_associations = pao1_similarity_scores.merge(
    pao1_expression_stats, left_index=True, right_index=True, how="left"
)
pa14_associations = pa14_similarity_scores.merge(
    pa14_expression_stats, left_index=True, right_index=True, how="left"
)

print(pao1_associations.shape)
pao1_associations.head()

print(pa14_associations.shape)
pa14_associations.head()

# Save
pao1_associations.to_csv(pao1_out_filename, sep="\t")
pa14_associations.to_csv(pa14_out_filename, sep="\t")

# ## Examine expression distribution
#
# One of the "most stable" core genes found were from the T6SS, which is surprising given that...
#
# We will plot the distribution of these genes to make sure that ____
#
# We manually selected these genes.

# ### Example of most stable core gene

# +
# tssC1 (T6SS) gene selected
pao1_most_id = "PA0084"
pa14_most_id = "PA14_01020"

# hcp1
pao1_most_id = "PA0085"
pa14_most_id = "PA14_01030"
# -

sns.displot(np.log10(pao1_expression.loc[pao1_most_id]))
sns.displot(np.log10(pa14_expression.loc[pa14_most_id]))

# ### Example of least stable core gene

# +
# Least stable core gene
# pao1_least_id = "PA3507"
# pa14_least_id = "PA14_10840"

# gloA2
pao1_least_id = "PA0710"
pa14_least_id = "PA14_55130"
# -

sns.displot(np.log10(pao1_expression.loc[pao1_least_id]))
sns.displot(np.log10(pa14_expression.loc[pa14_least_id]))
