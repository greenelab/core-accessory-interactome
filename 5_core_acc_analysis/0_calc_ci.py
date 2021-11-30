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

# # Calculate confidence interval
#
# This notebook generates the confidence interval for the plot in [stable_gene_relationships.ipynb](stable_gene_relationships.ipynb) notebook. Since this confidence interval is based on boostraping it takes a while to run so we wanted this in a separate notebook.
#
# Existing ci calculation either assumes normality or uses bootstrapping, but since we need to make adjustments to normalize our results we cannot use the out of the box bootstrapping
# https://stackoverflow.com/questions/46125182/is-seaborn-confidence-interval-computed-correctly

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import scipy
import pandas as pd
import numpy as np
import textwrap
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scripts import utils, paths, gene_relationships, annotations

random.seed(1)

# +
# User params
offset_to_bin = 10

use_operon = True
sum_increment_to_use = 1

n_boot = 50

ci = 0.95

# Output filename
pao1_figure_filename = (
    "PAO1_stablility_expression_relationships_operon_corrected_spell.svg"
)
pa14_figure_filename = (
    "PA14_stability_expression_relationships_operon_corrected_spell.svg"
)
# -

# ### Import gene ids
#

# +
# Import correlation matrix to get gene ids
pao1_corr_filename = paths.PAO1_CORR_LOG_SPELL
pa14_corr_filename = paths.PA14_CORR_LOG_SPELL

pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pa14_corr = pd.read_csv(pa14_corr_filename, sep="\t", index_col=0, header=0)
# -

# Make a dataframe with gene ids
pao1_membership = pd.DataFrame(data=[], index=pao1_corr.index)
print(pao1_membership.shape)
pao1_membership.head()

pa14_membership = pd.DataFrame(data=[], index=pa14_corr.index)
print(pa14_membership.shape)
pa14_membership.head()

# ### Import and format operon data

pao1_operon_filename = paths.PAO1_OPERON
pa14_operon_filename = paths.PA14_OPERON

pao1_operon = annotations.load_format_operons(pao1_operon_filename)
pa14_operon = annotations.load_format_operons(pa14_operon_filename)

if use_operon:
    pao1_operon_expression_to_use = pao1_operon
    pa14_operon_expression_to_use = pa14_operon
else:
    pao1_operon_expression_to_use = None
    pa14_operon_expression_to_use = None

# ### Map core/accessory labels to genes

# Read in expression data
pao1_expression_filename = paths.PAO1_COMPENDIUM
pa14_expression_filename = paths.PA14_COMPENDIUM

pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT

(
    pao1_arr,
    pa14_arr,
    pao1_core,
    pao1_acc,
    pa14_core,
    pa14_acc,
) = annotations.map_core_acc_annot(
    pao1_membership,
    pa14_membership,
    pao1_expression_filename,
    pa14_expression_filename,
    pao1_annot_filename,
    pa14_annot_filename,
)

# ## Find relationships using expression distance

# Correlation matrix files
pao1_corr_filename = paths.PAO1_CORR_LOG_SPELL
pa14_corr_filename = paths.PA14_CORR_LOG_SPELL

# Load correlation data
pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pa14_corr = pd.read_csv(pa14_corr_filename, sep="\t", index_col=0, header=0)

# +
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

# +
# Get most and least stable core genes
pao1_most_stable_genes = list(
    pao1_similarity_scores[pao1_similarity_scores["label"] == "most stable"].index
)
pao1_least_stable_genes = list(
    pao1_similarity_scores[pao1_similarity_scores["label"] == "least stable"].index
)

pa14_most_stable_genes = list(
    pa14_similarity_scores[pa14_similarity_scores["label"] == "most stable"].index
)
pa14_least_stable_genes = list(
    pa14_similarity_scores[pa14_similarity_scores["label"] == "least stable"].index
)
# -

# %%time
expression_dist_counts_pao1_most_ci = (
    gene_relationships.get_CI_expression_relationships(
        n_boot,
        pao1_corr,
        pao1_most_stable_genes,
        pao1_arr,
        offset_to_bin,
        pao1_operon_expression_to_use,
        sum_increment_to_use,
    )
)

# %%time
expression_dist_counts_pao1_least_ci = (
    gene_relationships.get_CI_expression_relationships(
        n_boot,
        pao1_corr,
        pao1_least_stable_genes,
        pao1_arr,
        offset_to_bin,
        pao1_operon_expression_to_use,
        sum_increment_to_use,
    )
)

# %%time
expression_dist_counts_pa14_most_ci = (
    gene_relationships.get_CI_expression_relationships(
        n_boot,
        pa14_corr,
        pa14_most_stable_genes,
        pa14_arr,
        offset_to_bin,
        pa14_operon_expression_to_use,
        sum_increment_to_use,
    )
)

# %%time
expression_dist_counts_pa14_least_ci = (
    gene_relationships.get_CI_expression_relationships(
        n_boot,
        pa14_corr,
        pa14_least_stable_genes,
        pa14_arr,
        offset_to_bin,
        pa14_operon_expression_to_use,
        sum_increment_to_use,
    )
)

expression_dist_counts_pao1_most_ci

expression_dist_counts_pao1_least_ci

# ## Calculate percentages
# Here we are taking the number of co-expressed core or accessory genes and normalizing by the number of most or least stable genes to get a percentage

# Get only columns with counts from bootstrapping
sampling_cols = [
    col for col in expression_dist_counts_pao1_most_ci.columns if "total" in col
]

# +
pao1_acc_most_ids = expression_dist_counts_pao1_most_ci.loc[
    expression_dist_counts_pao1_most_ci["gene type"] == "acc"
].index
pao1_core_most_ids = expression_dist_counts_pao1_most_ci.loc[
    expression_dist_counts_pao1_most_ci["gene type"] == "core"
].index

pao1_acc_least_ids = expression_dist_counts_pao1_least_ci.loc[
    expression_dist_counts_pao1_least_ci["gene type"] == "acc"
].index
pao1_core_least_ids = expression_dist_counts_pao1_least_ci.loc[
    expression_dist_counts_pao1_least_ci["gene type"] == "core"
].index

# +
pa14_acc_most_ids = expression_dist_counts_pa14_most_ci.loc[
    expression_dist_counts_pa14_most_ci["gene type"] == "acc"
].index
pa14_core_most_ids = expression_dist_counts_pa14_most_ci.loc[
    expression_dist_counts_pa14_most_ci["gene type"] == "core"
].index

pa14_acc_least_ids = expression_dist_counts_pa14_least_ci.loc[
    expression_dist_counts_pa14_least_ci["gene type"] == "acc"
].index
pa14_core_least_ids = expression_dist_counts_pa14_least_ci.loc[
    expression_dist_counts_pa14_least_ci["gene type"] == "core"
].index
# -

# Most stable PAO1
expression_dist_counts_pao1_most_ci.loc[pao1_acc_most_ids, sampling_cols] /= len(
    pao1_most_stable_genes
)
expression_dist_counts_pao1_most_ci.loc[pao1_core_most_ids, sampling_cols] /= len(
    pao1_most_stable_genes
)

# Least stable PAO1
expression_dist_counts_pao1_least_ci.loc[pao1_acc_least_ids, sampling_cols] /= len(
    pao1_least_stable_genes
)
expression_dist_counts_pao1_least_ci.loc[pao1_core_least_ids, sampling_cols] /= len(
    pao1_least_stable_genes
)

# Most stable PA14
expression_dist_counts_pa14_most_ci.loc[pa14_acc_most_ids, sampling_cols] /= len(
    pa14_most_stable_genes
)
expression_dist_counts_pa14_most_ci.loc[pa14_core_most_ids, sampling_cols] /= len(
    pa14_most_stable_genes
)

# Least stable PA14
expression_dist_counts_pa14_least_ci.loc[pa14_acc_least_ids, sampling_cols] /= len(
    pa14_least_stable_genes
)
expression_dist_counts_pa14_least_ci.loc[pa14_core_least_ids, sampling_cols] /= len(
    pa14_least_stable_genes
)

expression_dist_counts_pao1_most_ci

expression_dist_counts_pao1_least_ci

# ## Normalize by base percentage
# Here we want to normalize the percentage of co-expressed genes with the percent of accessory or core genes in the genome.

# Baseline/expected proportions for PAO1
pao1_total = len(pao1_core) + len(pao1_acc)
pao1_acc_expected = len(pao1_acc) / pao1_total
pao1_core_expected = len(pao1_core) / pao1_total
print("total pao1 genes", pao1_total)
print("pao1 acc baseline", pao1_acc_expected)
print("pao1 core baseline", pao1_core_expected)

# Baseline/expected proportions for PA14
pa14_total = len(pa14_core) + len(pa14_acc)
pa14_acc_expected = len(pa14_acc) / pa14_total
pa14_core_expected = len(pa14_core) / pa14_total
print("total pa14 genes", pa14_total)
print("pa14 acc baseline", pa14_acc_expected)
print("pa14 core baseline", pa14_core_expected)

# Most stable PAO1
expression_dist_counts_pao1_most_ci.loc[
    pao1_acc_most_ids, sampling_cols
] /= pao1_acc_expected
expression_dist_counts_pao1_most_ci.loc[
    pao1_core_most_ids, sampling_cols
] /= pao1_core_expected

# Least stable PAO1
expression_dist_counts_pao1_least_ci.loc[
    pao1_acc_least_ids, sampling_cols
] /= pao1_acc_expected
expression_dist_counts_pao1_least_ci.loc[
    pao1_core_least_ids, sampling_cols
] /= pao1_core_expected

# Most stable PA14
expression_dist_counts_pa14_most_ci.loc[
    pa14_acc_most_ids, sampling_cols
] /= pa14_acc_expected
expression_dist_counts_pa14_most_ci.loc[
    pa14_core_most_ids, sampling_cols
] /= pa14_core_expected

# Least stable PA14
expression_dist_counts_pa14_least_ci.loc[
    pa14_acc_least_ids, sampling_cols
] /= pa14_acc_expected
expression_dist_counts_pa14_least_ci.loc[
    pa14_core_least_ids, sampling_cols
] /= pa14_core_expected

expression_dist_counts_pao1_most_ci

expression_dist_counts_pao1_least_ci

# ## Get quantiles

alpha = 1 - ci
lower = alpha / 2
upper = 1 - (alpha / 2)

pao1_most_ci_ranges = expression_dist_counts_pao1_most_ci.quantile(
    [lower, upper], axis=1
)
pao1_least_ci_ranges = expression_dist_counts_pao1_least_ci.quantile(
    [lower, upper], axis=1
)

pa14_most_ci_ranges = expression_dist_counts_pa14_most_ci.quantile(
    [lower, upper], axis=1
)
pa14_least_ci_ranges = expression_dist_counts_pa14_least_ci.quantile(
    [lower, upper], axis=1
)

# ## Format
#
# Merge with starting df with corr

# +
pao1_most_ci = expression_dist_counts_pao1_most_ci.merge(
    pao1_most_ci_ranges.T, left_index=True, right_index=True
).drop(sampling_cols, axis=1)
pao1_least_ci = expression_dist_counts_pao1_least_ci.merge(
    pao1_least_ci_ranges.T, left_index=True, right_index=True
).drop(sampling_cols, axis=1)

pao1_most_ci.columns = ["offset", "gene type", "ymin", "ymax"]
pao1_least_ci.columns = ["offset", "gene type", "ymin", "ymax"]

# +
pa14_most_ci = expression_dist_counts_pa14_most_ci.merge(
    pa14_most_ci_ranges.T, left_index=True, right_index=True
).drop(sampling_cols, axis=1)
pa14_least_ci = expression_dist_counts_pa14_least_ci.merge(
    pa14_least_ci_ranges.T, left_index=True, right_index=True
).drop(sampling_cols, axis=1)

pa14_most_ci.columns = ["offset", "gene type", "ymin", "ymax"]
pa14_least_ci.columns = ["offset", "gene type", "ymin", "ymax"]
# -

pao1_most_ci.head()

pa14_most_ci.head()

# +
# Save
pao1_most_ci.to_csv("pao1_most_ci.tsv", sep="\t")
pao1_least_ci.to_csv("pao1_least_ci.tsv", sep="\t")

pa14_most_ci.to_csv("pa14_most_ci.tsv", sep="\t")
pa14_least_ci.to_csv("pa14_least_ci.tsv", sep="\t")
