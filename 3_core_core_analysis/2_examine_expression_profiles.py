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
pao1_out_filename = "pao1_core_similarity_expression_stats_spell.tsv"
pa14_out_filename = "pa14_core_similarity_expression_stats_spell.tsv"

# +
# Load transcriptional similarity df
pao1_similarity_scores_filename = "pao1_similarity_scores_spell.tsv"
pa14_similarity_scores_filename = "pa14_similarity_scores_spell.tsv"

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

# Expression matrices are sample x gene
pao1_expression = pd.read_csv(pao1_expression_filename, sep="\t", header=0, index_col=0)
pa14_expression = pd.read_csv(pa14_expression_filename, sep="\t", header=0, index_col=0)

pao1_expression.head()

pao1_expression.tail()

# Get distribution statistics
pao1_expression_stats = pao1_expression.describe()
pa14_expression_stats = pa14_expression.describe()

pao1_expression_stats.head()

pa14_expression_stats.head()

# Format statistic data
pao1_expression_stats = pao1_expression_stats.T
pao1_expression_stats = pao1_expression_stats.drop(columns=["count"])
pao1_expression_stats = pao1_expression_stats.rename(
    columns={
        "mean": "mean expression",
        "std": "standard deviation expression",
        "25%": "25% expression",
        "50%": "50% expression",
        "75%": "75% expression",
        "min": "min expression",
        "max": "max expression",
    }
)
pao1_expression_stats["variance expression"] = (
    pao1_expression_stats["standard deviation expression"] ** 2
)
pao1_expression_stats["range expression"] = (
    pao1_expression_stats["max expression"] - pao1_expression_stats["min expression"]
)

pa14_expression_stats = pa14_expression_stats.T
pa14_expression_stats = pa14_expression_stats.drop(columns=["count"])
pa14_expression_stats = pa14_expression_stats.rename(
    columns={
        "mean": "mean expression",
        "std": "standard deviation expression",
        "25%": "25% expression",
        "50%": "50% expression",
        "75%": "75% expression",
        "min": "min expression",
        "max": "max expression",
    }
)
pa14_expression_stats["variance expression"] = (
    pa14_expression_stats["standard deviation expression"] ** 2
)
pa14_expression_stats["range expression"] = (
    pa14_expression_stats["max expression"] - pa14_expression_stats["min expression"]
)

pao1_expression_stats.head()

pa14_expression_stats.head()

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
# One of the "most stable" core genes found were from the T6SS, which is surprising given this pathway allows for inter-strain warfare and so we’d expect genes within this pathway to vary across strains.
#
# We will plot the distribution of these genes to make sure that the reason these T6SS genes are found to be stable is because all the genes are "off". Based on the plots below, this doesn't look to be the explanation for why T6SS genes are found to be stable across strains. Genes are expressed (> 1.0 log10 expression = 10 normalized counts)
#
# We manually selected these genes.

# ### Example of most stable core gene

# +
# tssC1 (T6SS) gene selected
# pao1_most_id = "PA0084"
# pa14_most_id = "PA14_01020"

# hcp1 (T6SS)
# pao1_most_id = "PA0085"
# pa14_most_id = "PA14_01030"

# tssF1 (T6SS)
# pao1_most_id = "PA0088"
# pa14_most_id = "PA14_01070"

# pscC (T3SS)
# pao1_most_id = "PA1716"
# pa14_most_id = "PA14_42350"

# pscF (T3SS)
pao1_most_id = "PA1719"
pa14_most_id = "PA14_42310"
# -

sns.displot(np.log10(pao1_expression[pao1_most_id] + 1))
sns.displot(np.log10(pa14_expression[pa14_most_id] + 1))

# ### Example of least stable core gene

# +
# Least stable core gene
# pao1_least_id = "PA4685"
# pa14_least_id = "PA14_61980"

pao1_least_id = "PA2458"
pa14_least_id = "PA14_32830"
# -

sns.displot(np.log10(pao1_expression[pao1_least_id] + 1))
sns.displot(np.log10(pa14_expression[pa14_least_id] + 1))
