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

# # Relationships using expression distance
#
# This notebook is performing the same analysis as seen in [all_gene_relationships.ipynb](archive/all_gene_relationships.ipynb), where we are examining who is related to who. Previously we started with an accessory gene and asked: is the highest correlated gene another accessory gene or a core gene? For this analysis, we are starting with the most stable core genes and asking the same question: is the highest correlated gene core or accessory?
#
# Note: We do not have the genome location metric here because this would require a significant effort to figure out how to modify the existing code to only focus on a subset of genes.

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

# Output filename
pao1_figure_filename = (
    "PAO1_stablility_expression_relationships_operon_corrected_spell.svg"
)
pa14_figure_filename = (
    "PA14_stability_expression_relationships_operon_corrected_spell.svg"
)
# -

# ### Import gene ids

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

print(pao1_operon.shape)
pao1_operon.head()

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

print(pao1_arr.shape)
pao1_arr.head()

pao1_arr.tail()

print(pa14_arr.shape)
pa14_arr.head()

pa14_arr.tail()

# +
# Fill in index of operon_df to include all genes
all_pao1_gene_ids = pao1_arr.index
all_pa14_gene_ids = pa14_arr.index

# Get missing gene ids
missing_pao1_gene_ids = set(all_pao1_gene_ids).difference(pao1_operon.index)
missing_pa14_gene_ids = set(all_pa14_gene_ids).difference(pa14_operon.index)

# Make dataframe with missing gene ids with np.nan values for operon_name
missing_pao1_gene_df = pd.DataFrame(
    data=np.nan, index=list(missing_pao1_gene_ids), columns=["operon_name"]
)
missing_pa14_gene_df = pd.DataFrame(
    data=np.nan, index=list(missing_pa14_gene_ids), columns=["operon_name"]
)

pao1_operon_genome_dist = pao1_operon.append(missing_pao1_gene_df)
pa14_operon_genome_dist = pa14_operon.append(missing_pa14_gene_df)

pao1_operon_genome_dist = pao1_operon_genome_dist.loc[all_pao1_gene_ids]
pa14_operon_genome_dist = pa14_operon_genome_dist.loc[all_pa14_gene_ids]
# -

print(pao1_operon_genome_dist.shape)
pao1_operon_genome_dist.tail()

print(pa14_operon_genome_dist.shape)
pa14_operon_genome_dist.tail()

if use_operon:
    pao1_operon_genome_to_use = pao1_operon_genome_dist
    pa14_operon_genome_to_use = pa14_operon_genome_dist
else:
    pao1_operon_genome_to_use = None
    pa14_operon_genome_to_use = None

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
expression_dist_counts_pao1_most = (
    gene_relationships.get_relationship_in_expression_space(
        pao1_corr,
        pao1_most_stable_genes,
        pao1_arr,
        offset_to_bin,
        pao1_operon_expression_to_use,
        sum_increment_to_use,
    )
)

# +
# %%time
expression_dist_counts_pao1_most_ci = (
    gene_relationships.get_CI_expression_relationships(
        5,
        pao1_corr,
        pao1_most_stable_genes,
        pao1_arr,
        offset_to_bin,
        pao1_operon_expression_to_use,
        sum_increment_to_use,
    )
)

# TO DO:
# To finish in the next PR
# Check existing method to calculate ci in sns.barplot
# https://stackoverflow.com/questions/46125182/is-seaborn-confidence-interval-computed-correctly
# Customize ci with https://stackoverflow.com/questions/52767423/is-it-possible-to-input-values-for-confidence-interval-error-bars-on-seaborn-ba
# Add this for other relationships
# -

# %%time
expression_dist_counts_pao1_least = (
    gene_relationships.get_relationship_in_expression_space(
        pao1_corr,
        pao1_least_stable_genes,
        pao1_arr,
        offset_to_bin,
        pao1_operon_expression_to_use,
        sum_increment_to_use,
    )
)

# %%time
expression_dist_counts_pa14_most = (
    gene_relationships.get_relationship_in_expression_space(
        pa14_corr,
        pa14_most_stable_genes,
        pa14_arr,
        offset_to_bin,
        pa14_operon_expression_to_use,
        sum_increment_to_use,
    )
)

# %%time
expression_dist_counts_pa14_least = (
    gene_relationships.get_relationship_in_expression_space(
        pa14_corr,
        pa14_least_stable_genes,
        pa14_arr,
        offset_to_bin,
        pa14_operon_expression_to_use,
        sum_increment_to_use,
    )
)

expression_dist_counts_pao1_most.head()

expression_dist_counts_pao1_least.head()

expression_dist_counts_pa14_most.head()

expression_dist_counts_pa14_least.head()

# ### Format data for plotting
#
# Here we will calculate the proportion of gene types per offset and then normalize by the proportion of core and accessory genes. This will return an oddsratio type value where if the value is >1 than the proportion of genes of that type are more than expected.

# Calculate the percentages per offset
expression_dist_counts_pao1_most["percent"] = expression_dist_counts_pao1_most[
    "total"
] / len(pao1_most_stable_genes)
expression_dist_counts_pao1_least["percent"] = expression_dist_counts_pao1_least[
    "total"
] / len(pao1_least_stable_genes)

expression_dist_counts_pa14_most["percent"] = expression_dist_counts_pa14_most[
    "total"
] / len(pa14_most_stable_genes)
expression_dist_counts_pa14_least["percent"] = expression_dist_counts_pa14_least[
    "total"
] / len(pa14_least_stable_genes)

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

# +
# Normalize by baseline PAO1 most stable
pao1_acc_most_ids = expression_dist_counts_pao1_most.loc[
    expression_dist_counts_pao1_most["gene type"] == "acc"
].index
pao1_core_most_ids = expression_dist_counts_pao1_most.loc[
    expression_dist_counts_pao1_most["gene type"] == "core"
].index

expression_dist_counts_pao1_most.loc[pao1_acc_most_ids, "normalized"] = (
    expression_dist_counts_pao1_most.loc[pao1_acc_most_ids, "percent"]
    / pao1_acc_expected
)
expression_dist_counts_pao1_most.loc[pao1_core_most_ids, "normalized"] = (
    expression_dist_counts_pao1_most.loc[pao1_core_most_ids, "percent"]
    / pao1_core_expected
)

# CI
# expression_dist_counts_pao1_most_ci[
#    ["total_0", "total_1", "total_2", "total_3", "total_4"]
# ] /= pao1_acc_expected
# Get 95% range
# pao1_most_ci_ranges = expression_dist_counts_pao1_most_ci.quantile(
#    [0.025, 0.975], axis=1
# )
# -

expression_dist_counts_pao1_most_ci

# Get percent
expression_dist_counts_pao1_most_ci.loc[
    pao1_acc_most_ids, ["total_0", "total_1", "total_2", "total_3", "total_4"]
] /= len(pao1_most_stable_genes)

expression_dist_counts_pao1_most_ci

# Normalize
expression_dist_counts_pao1_most_ci.loc[
    pao1_acc_most_ids, ["total_0", "total_1", "total_2", "total_3", "total_4"]
] /= pao1_acc_expected

expression_dist_counts_pao1_most_ci

# Get 95% range
pao1_most_ci_ranges = expression_dist_counts_pao1_most_ci.quantile(
    [0.025, 0.975], axis=1
)

pao1_most_ci_ranges

# +
## TESTING ABOVE

# +
# Normalize by baseline PAO1 least stable
pao1_acc_least_ids = expression_dist_counts_pao1_least.loc[
    expression_dist_counts_pao1_least["gene type"] == "acc"
].index
pao1_core_least_ids = expression_dist_counts_pao1_least.loc[
    expression_dist_counts_pao1_least["gene type"] == "core"
].index

expression_dist_counts_pao1_least.loc[pao1_acc_least_ids, "normalized"] = (
    expression_dist_counts_pao1_least.loc[pao1_acc_least_ids, "percent"]
    / pao1_acc_expected
)
expression_dist_counts_pao1_least.loc[pao1_core_least_ids, "normalized"] = (
    expression_dist_counts_pao1_least.loc[pao1_core_least_ids, "percent"]
    / pao1_core_expected
)

# +
# Normalize by baseline PA14 most stable
pa14_acc_most_ids = expression_dist_counts_pa14_most.loc[
    expression_dist_counts_pa14_most["gene type"] == "acc"
].index
pa14_core_most_ids = expression_dist_counts_pao1_most.loc[
    expression_dist_counts_pa14_most["gene type"] == "core"
].index

expression_dist_counts_pa14_most.loc[pa14_acc_most_ids, "normalized"] = (
    expression_dist_counts_pa14_most.loc[pa14_acc_most_ids, "percent"]
    / pa14_acc_expected
)
expression_dist_counts_pa14_most.loc[pa14_core_most_ids, "normalized"] = (
    expression_dist_counts_pa14_most.loc[pa14_core_most_ids, "percent"]
    / pa14_core_expected
)

# +
# Normalize by baseline PA14 least stable
pa14_acc_least_ids = expression_dist_counts_pa14_least.loc[
    expression_dist_counts_pa14_least["gene type"] == "acc"
].index
pa14_core_least_ids = expression_dist_counts_pa14_least.loc[
    expression_dist_counts_pa14_least["gene type"] == "core"
].index

expression_dist_counts_pa14_least.loc[pa14_acc_least_ids, "normalized"] = (
    expression_dist_counts_pa14_least.loc[pa14_acc_least_ids, "percent"]
    / pa14_acc_expected
)
expression_dist_counts_pa14_least.loc[pa14_core_least_ids, "normalized"] = (
    expression_dist_counts_pa14_least.loc[pa14_core_least_ids, "percent"]
    / pa14_core_expected
)

# +
# Combine PAO1 dataframes
expression_dist_counts_pao1_most.loc[pao1_acc_most_ids, "label"] = "most stable acc"
expression_dist_counts_pao1_most.loc[pao1_core_most_ids, "label"] = "most stable core"
expression_dist_counts_pao1_least.loc[pao1_acc_least_ids, "label"] = "least stable acc"
expression_dist_counts_pao1_least.loc[
    pao1_core_least_ids, "label"
] = "least stable core"

expression_dist_counts_pao1_all = pd.concat(
    [expression_dist_counts_pao1_most, expression_dist_counts_pao1_least]
)
# -

expression_dist_counts_pao1_all

# +
# Combine PA14 dataframes
expression_dist_counts_pa14_most.loc[pa14_acc_most_ids, "label"] = "most stable acc"
expression_dist_counts_pa14_most.loc[pa14_core_most_ids, "label"] = "most stable core"
expression_dist_counts_pa14_least.loc[pa14_acc_least_ids, "label"] = "least stable acc"
expression_dist_counts_pa14_least.loc[
    pa14_core_least_ids, "label"
] = "least stable core"

expression_dist_counts_pa14_all = pd.concat(
    [expression_dist_counts_pa14_most, expression_dist_counts_pa14_least]
)
# -

expression_dist_counts_pa14_all

# ### Plot

# +
# Plot PAO1 trends
plt.figure(figsize=(10, 8))

pao1_subset = expression_dist_counts_pao1_all[
    (expression_dist_counts_pao1_all["gene type"] == "acc")
    & (expression_dist_counts_pao1_all["label"] == "most stable acc")
]

fig = sns.barplot(
    data=pao1_subset,
    x="offset",
    y="normalized",
    hue="label",
    hue_order=[
        "most stable acc",
        "least stable acc",
        # "most stable core",
        # "least stable core",
    ],
    palette={
        "most stable acc": "#21368B",
        "least stable acc": "#A6AED0",
        # "most stable core": "#F8744C",
        # "least stable core": "#FCC7B7",
    },
    estimator=np.mean,
    ci=90,
    n_boot=5,
    errcolor="red",
    errwidth=2,
)
plt.errorbar(
    x=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "10+"],
    y=pao1_subset["normalized"],
    yerr=(
        (
            pao1_most_ci_ranges.loc[0.975, pao1_acc_most_ids]
            - pao1_most_ci_ranges.loc[0.025, pao1_acc_most_ids]
        )
    ),
)

plt.axhline(y=1.0, color="black", linestyle="--")
fig.set_title("Who are most/least stable core genes related to (PAO1)", fontsize=16)
fig.set_ylabel(
    r"Fold change" + "\n" + "(% acc genes co-express/% acc genes in genome)",
    fontsize=14,
)
fig.set_xlabel("Rank correlation in expression space", fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 0.6), loc=2, borderaxespad=0.0, fontsize=12)

# new_labels = ["Most stable gene", "Least stable gene"]
# for t, l in zip(fig._legend.texts, new_labels): t.set_text(l)
# -

"""pao1_subset = expression_dist_counts_pao1_all[
    (expression_dist_counts_pao1_all["gene type"] == "acc") &
    (expression_dist_counts_pao1_all["label"] == "most stable acc")
]
pao1_subset["normalized"]"""

# +
# Plot PA14 trends
plt.figure(figsize=(10, 8))

fig2 = sns.barplot(
    data=expression_dist_counts_pa14_all[
        expression_dist_counts_pa14_all["gene type"] == "acc"
    ],
    x="offset",
    y="normalized",
    hue="label",
    hue_order=[
        "most stable acc",
        "least stable acc",
        # "most stable core",
        # "least stable core",
    ],
    palette={
        "most stable acc": "#21368B",
        "least stable acc": "#A6AED0",
        # "most stable core": "#F8744C",
        # "least stable core": "#FCC7B7",
    },
)
plt.axhline(y=1.0, color="black", linestyle="--")
fig2.set_title("Who are most/least stable core genes related to (PA14)", fontsize=16)
fig2.set_ylabel(
    r"Fold change" + "\n" + "(% acc genes co-express/% acc genes in genome)",
    fontsize=14,
)
fig2.set_xlabel("Rank correlation in expression space", fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 0.6), loc=2, borderaxespad=0.0, fontsize=12)
# -

"""# Save figures using operons
# Save figures not using operons
# Save figure with rolling sum and operons
# Save figure with rolling sum not using operons
fig.figure.savefig(
    pao1_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

fig2.figure.savefig(
    pa14_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)"""

# +
# TO DO:
# Legend label
# Fix CI addition
# plt.errorbar(
#    x=['1','2','3','4','5','6','7','8','9','10','10+'],
#    yerr=pao1_most_ci_ranges
# )
# -

# **Takeaway:**
#
# * Least stable core genes have more accessory gene neighbors compared to most stable core genes
# * Previous evidence found that insertion sequences (type of accessory gene) can change the expression of existing genes once it is integrated into the genome. So perhaps these least stable core genes transcriptional behavior is modified by the accessory genes.
