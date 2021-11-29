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
from plotnine import (
    ggplot,
    labs,
    geom_hline,
    geom_bar,
    geom_errorbar,
    positions,
    aes,
    ggsave,
    theme_bw,
    theme,
    facet_wrap,
    scale_fill_manual,
    scale_x_discrete,
    xlim,
    guides,
    guide_legend,
    element_blank,
    element_text,
    element_rect,
    element_line,
    coords,
)

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
# -

# Combine PAO1 dataframes
expression_dist_counts_pao1_most.loc[pao1_acc_most_ids, "label"] = "most stable acc"
expression_dist_counts_pao1_most.loc[pao1_core_most_ids, "label"] = "most stable core"
expression_dist_counts_pao1_least.loc[pao1_acc_least_ids, "label"] = "least stable acc"
expression_dist_counts_pao1_least.loc[
    pao1_core_least_ids, "label"
] = "least stable core"

# Combine PA14 dataframes
expression_dist_counts_pa14_most.loc[pa14_acc_most_ids, "label"] = "most stable acc"
expression_dist_counts_pa14_most.loc[pa14_core_most_ids, "label"] = "most stable core"
expression_dist_counts_pa14_least.loc[pa14_acc_least_ids, "label"] = "least stable acc"
expression_dist_counts_pa14_least.loc[
    pa14_core_least_ids, "label"
] = "least stable core"


# ### Add confidence interval

# +
# Import confidence interval data
pao1_most_ci = pd.read_csv("pao1_most_ci.tsv", sep="\t", index_col=0, header=0)
pao1_least_ci = pd.read_csv("pao1_least_ci.tsv", sep="\t", index_col=0, header=0)

pa14_most_ci = pd.read_csv("pa14_most_ci.tsv", sep="\t", index_col=0, header=0)
pa14_least_ci = pd.read_csv("pa14_least_ci.tsv", sep="\t", index_col=0, header=0)
# -

expression_dist_counts_pao1_most = expression_dist_counts_pao1_most.merge(
    pao1_most_ci[["ymin", "ymax"]], left_index=True, right_index=True
)
expression_dist_counts_pao1_least = expression_dist_counts_pao1_least.merge(
    pao1_least_ci[["ymin", "ymax"]], left_index=True, right_index=True
)

expression_dist_counts_pa14_most = expression_dist_counts_pa14_most.merge(
    pa14_most_ci[["ymin", "ymax"]], left_index=True, right_index=True
)
expression_dist_counts_pa14_least = expression_dist_counts_pa14_least.merge(
    pa14_least_ci[["ymin", "ymax"]], left_index=True, right_index=True
)

expression_dist_counts_pao1_all = pd.concat(
    [expression_dist_counts_pao1_most, expression_dist_counts_pao1_least]
)

expression_dist_counts_pao1_all

expression_dist_counts_pa14_all = pd.concat(
    [expression_dist_counts_pa14_most, expression_dist_counts_pa14_least]
)

expression_dist_counts_pa14_all

# ### Plot

# +
pao1_subset = expression_dist_counts_pao1_all[
    (expression_dist_counts_pao1_all["gene type"] == "acc")
]
pao1_subset["offset"] = list(pao1_subset["offset"].astype("str"))
x_ticks = ["+10", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]

fig_pao1 = (
    ggplot(pao1_subset, aes(x="offset", y="normalized", fill="label"))
    + geom_bar(stat="identity", position="dodge", width=0.8)
    + geom_errorbar(
        pao1_subset,
        aes(x="offset", ymin="ymin", ymax="ymax"),
        position=positions.position_dodge(0.8),
        color="darkgrey",
    )
    + geom_hline(aes(yintercept=1.0), linetype="dashed")
    + labs(
        x="Rank co-expression",
        y=r"Fold change" + "\n" + "(% acc genes co-express/% acc genes in genome)",
        title="Who are most/least stable core genes related to (PAO1)",
    )
    + theme(
        plot_background=element_rect(fill="white"),
        panel_background=element_rect(fill="white"),
        panel_grid_major_x=element_line(color="lightgrey"),
        panel_grid_major_y=element_line(color="lightgrey"),
        axis_line=element_line(color="grey"),
        legend_key=element_rect(fill="white", colour="white"),
        legend_title=element_blank(),
        legend_text=element_text(family="sans-serif", size=12),
        plot_title=element_text(family="sans-serif", size=15),
        axis_text=element_text(family="sans-serif", size=12),
        axis_title=element_text(family="sans-serif", size=10),
    )
    + scale_fill_manual(
        values=["#21368B", "#A6AED0"],
        labels=[
            "Accessory genes related to least stable core genes",
            "Accessory genes related to most stable core genes",
        ],
    )
    + scale_x_discrete(limits=x_ticks, labels=x_ticks)
)
print(fig_pao1)

# +
pa14_subset = expression_dist_counts_pa14_all[
    (expression_dist_counts_pa14_all["gene type"] == "acc")
]
pa14_subset["offset"] = list(pa14_subset["offset"].astype("str"))

fig_pa14 = (
    ggplot(pa14_subset, aes(x="offset", y="normalized", fill="label"))
    + geom_bar(stat="identity", position="dodge", width=0.8)
    + geom_errorbar(
        pa14_subset,
        aes(x="offset", ymin="ymin", ymax="ymax"),
        position=positions.position_dodge(0.8),
        color="darkgrey",
    )
    + geom_hline(aes(yintercept=1.0), linetype="dashed")
    + labs(
        x="Rank co-expression",
        y=r"Fold change" + "\n" + "(% acc genes co-express/% acc genes in genome)",
        title="Who are most/least stable core genes related to (PA14)",
    )
    + theme(
        plot_background=element_rect(fill="white"),
        panel_background=element_rect(fill="white"),
        panel_grid_major_x=element_line(color="lightgrey"),
        panel_grid_major_y=element_line(color="lightgrey"),
        axis_line=element_line(color="grey"),
        legend_key=element_rect(fill="white", colour="white"),
        legend_title=element_blank(),
        legend_text=element_text(family="sans-serif", size=12),
        plot_title=element_text(family="sans-serif", size=15),
        axis_text=element_text(family="sans-serif", size=12),
        axis_title=element_text(family="sans-serif", size=10),
    )
    + scale_fill_manual(
        values=["#21368B", "#A6AED0"],
        labels=[
            "Accessory genes related to least stable core genes",
            "Accessory genes related to most stable core genes",
        ],
    )
    + scale_x_discrete(limits=x_ticks, labels=x_ticks)
)
print(fig_pa14)

# +
# Calculate statistical test between the distribution of the top 10 co-expressed
# genes related to the least stable vs the most stable core genes
# Test: mean number of co-expressed accessory genes in least stable group vs mean number of
# co-expressed accessory genes in most stable group
# (compare dark blue and light blue bars)

pao1_least_df = pao1_subset[pao1_subset["label"] == "least stable acc"]
pao1_least_df = pao1_least_df[pao1_least_df.offset != "+10"]
pao1_least_vals = pao1_least_df["normalized"].values

pao1_most_df = pao1_subset[pao1_subset["label"] == "most stable acc"]
pao1_most_df = pao1_most_df[pao1_most_df.offset != "+10"]
pao1_most_vals = pao1_most_df["normalized"].values

# Independent t-test
# Test the null hypothesis such that the means of two populations are equal
(pao1_stats, pao1_pvalue) = scipy.stats.ttest_ind(pao1_least_vals, pao1_most_vals)
print(pao1_stats, pao1_pvalue)

# Non-parametric test
# nonparametric test of the null hypothesis that, for randomly selected values X and Y from two populations,
# the probability of X being greater than Y is equal to the probability of Y being greater than X.
(pao1_stats, pao1_pvalue) = scipy.stats.mannwhitneyu(pao1_least_vals, pao1_most_vals)
print(pao1_stats, pao1_pvalue)

# +
pa14_least_df = pa14_subset[pa14_subset["label"] == "least stable acc"]
pa14_least_df = pa14_least_df[pa14_least_df.offset != "+10"]
pa14_least_vals = pa14_least_df["normalized"].values

pa14_most_df = pa14_subset[pa14_subset["label"] == "most stable acc"]
pa14_most_df = pa14_most_df[pa14_most_df.offset != "+10"]
pa14_most_vals = pa14_most_df["normalized"].values

# Independent t-test
(pa14_stats, pa14_pvalue) = scipy.stats.ttest_ind(pa14_least_vals, pa14_most_vals)
print(pa14_stats, pa14_pvalue)

# Non-parametric test
(pa14_stats, pa14_pvalue) = scipy.stats.mannwhitneyu(pa14_least_vals, pa14_most_vals)
print(pa14_stats, pa14_pvalue)
# -

# Based on the bar plots we can be confident in our trend (as seen by the confidence intervals) that least stable genes are more co-expressed with accessory genes compared to most stable genes. This difference between least and most stable genes is further quantified by the t-test comparing the distribution of accessory genes related least vs most genes.

ggsave(plot=fig_pao1, filename=pao1_figure_filename, device="svg", dpi=300)
ggsave(plot=fig_pa14, filename=pa14_figure_filename, device="svg", dpi=300)

# **Takeaway:**
#
# * Least stable core genes have more accessory gene neighbors compared to most stable core genes
# * Previous evidence found that insertion sequences (type of accessory gene) can change the expression of existing genes once it is integrated into the genome. So perhaps these least stable core genes transcriptional behavior is modified by the accessory genes.
