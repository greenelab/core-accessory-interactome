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

# # Compare SPELL vs counts correlation
#
# This notebook performs an experiment to determine which correlation matrix we should use
# 1. Correlation of the MR counts expression matrix
# 2. Correlation of the MR counts that have been processed using SPELL
#
# The correlation of the counts matrix relates how similar a pair of genes are based on their expression profiles - **relates genes over samples**.
# High correlation means that a pair of genes have a similar expression profiles - i.e. similar levels of expression across samples/contexts, so genes both have low expression in the same samples and high expression in the same samples.
# * Pro: Easy to interpret
# * Con: Many gene pairs found to have a high correlation because many genes are related to the same pathway have the similar expression profiles. This is consistent with [Myers et al.](https://link.springer.com/article/10.1186/1471-2164-7-187), who found that there can be an over-representation of genes associated with the same pathway (i.e. a large fraction of gene pairs represent ribosomal relationships). This very prominent signal makes it difficult to detect other signals. Figure 1C demonstrates that a large fraction of gene pairs are ribosomal relationships - in the top 0.1% most co-expressed genes, 99% belong to the ribosome pathway. Furthermore, protein function prediction based on co-expression drop dramatically after removing the ribisome pathway (Figure 1A, B).
#
# To try to remove this very dominant global signal in the data. Here we are applying dimensionality reduction techniques in addition to scaling the data using a method called SPELL.
# The correlation of the SPELL matrix relates genes based on the gene coefficient matrix - **relate genes over their contribution to singular vectors (linear combination of genes - linear relationship between genes)**.
# High correlation means that a pair of genes contributes similarly to a singular vector, which are the axes pointing in the direction of the spread of the data and capture how genes are related to each other
# * Pro: Gene contributions are more balanced so that redundant signals (i.e. many genes from the same pathway - genes that vary together) are represented by a few SVs as opposed to many samples. More balanced also means that more subtle signals can be amplified (i.e. genes related by a smaller pathway are also captured by a few SVs)
# * Con: Can amplify noise - i.e. an SV that corresponds to some technical source of variability now has a similar weight to other real signals
#
# For more information comparing using counts vs SPELL-processing see: https://docs.google.com/presentation/d/18E0boNODJaxP-YYNIlccrh0kASbc7bapQBMovOX62jw/edit#slide=id.gf9d09c6be6_0_0

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import paths

# ## User params

# +
# Threshold to use to define edges between genes
# Top X% of genes are used
top_percent = 0.01

# Regulon file
regulon_filename = "gene_sets_refs.csv"
# -

# ## Load correlation matrix

pao1_corr_filename = paths.PAO1_CORR_RAW
pao1_corr_spell_filename = paths.PAO1_CORR_LOG_SPELL

# Load correlation data
pao1_corr_counts = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pao1_corr_spell = pd.read_csv(pao1_corr_spell_filename, sep="\t", index_col=0, header=0)


# ## Make edge matrix
#
# Convert correlation matrix of continuous values to an adjacency matrix with 1's if the correlation between a pair of genes exceeds the user defined threshold and therefore indicates if an edge exits those pair of genes.

# Get threshold to use based on percentage
def get_corr_threshold(corr_df, top_percent):
    # Since we are using the distribution of scores to determine the threshold
    # we need to remove duplicates and also the diagonal values
    # Here we get lower triangular matrix values only
    tril_corr_df = corr_df.where(~np.triu(np.ones(corr_df.shape)).astype(np.bool))

    # Flatten dataframe
    flat_corr_df = tril_corr_df.stack().reset_index()
    flat_corr_df.columns = ["gene_1", "gene_2", "corr_value"]

    # Get quantile
    # TO DO:Take abs????
    threshold = flat_corr_df.quantile(1 - top_percent)["corr_value"]
    print("correlation threshold: ", threshold)

    # Verify that number of gene pairs above the threshold
    # is approximately equal to the `top_percent`
    total_genes = flat_corr_df.shape[0]
    num_genes_above = flat_corr_df[flat_corr_df["corr_value"] > threshold].shape[0]
    percent_genes_above = num_genes_above / total_genes
    print("percent of pairs exceeding threshold: ", percent_genes_above)

    return threshold


pao1_corr_counts_threshold = get_corr_threshold(pao1_corr_counts, top_percent)

pao1_corr_spell_threshold = get_corr_threshold(pao1_corr_spell, top_percent)

# Create adjacency matrix using threshold defined above
# The adjacency matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by an edge
# TO DO:abs?????
pao1_counts_adj = (pao1_corr_counts > pao1_corr_counts_threshold).astype(float)
pao1_spell_adj = (pao1_corr_spell > pao1_corr_spell_threshold).astype(float)

pao1_counts_adj.head()

pao1_spell_adj.head()

# ## Evaluate relationships captured
#
# Given a set of known gene regulons, that we expect to be clustered together, let's calculate we want to evaluate which correlation matrix captures these geneset relationships better. To determine this we calculate the percentage of within (i.e. edges connecting two genes within the geneset) compared to the percentage of edges between(i.e.edges connecting a gene within the geneset and some other external gene).

# Load regulon data
regulon_df = pd.read_csv(regulon_filename, header=0, index_col=0)

# Format regulon data
regulon_df["Genes"] = regulon_df["Genes"].str.split(";").apply(list)

print(regulon_df.shape)
regulon_df.head()


def compare_within_btwn_edge(adj_df, regulon_df):
    """
    For each input regulon/geneset, this function calculates the
    percentage of edges within genes in the regulon vs the
    percentage of edges between genes in the regulon and some other
    gene.

    This function returns a dataframe containing the percentages for
    each regulon
    """

    # Loop through each regulon
    rows = []
    for regulon_name in regulon_df.index:
        geneset = regulon_df.loc[regulon_name, "Genes"]
        print(geneset)
        print(len(geneset))

        # Since there are some gene ids that are from PA14
        # We will take the intersection
        geneset_processed = set(geneset).intersection(adj_df.index)

        # Get within edges
        within_df = adj_df.loc[geneset_processed, geneset_processed]
        tril_within_df = within_df.where(
            ~np.triu(np.ones(within_df.shape)).astype(np.bool)
        )
        # print(tril_within_df)

        flat_within_df = tril_within_df.stack().reset_index()
        flat_within_df.columns = ["gene_1", "gene_2", "edge"]
        total_within_pairs = flat_within_df.shape[0]
        # print(flat_within_df)

        # count the number of within edges
        num_within_edges = flat_within_df["edge"].sum()
        print(num_within_edges)

        # Get between edges
        not_geneset = set(adj_df.index).difference(geneset_processed)
        between_df = adj_df.loc[geneset_processed, not_geneset]
        # tril_between_df = between_df.where(~np.triu(np.ones(between_df.shape)).astype(np.bool))

        flat_between_df = between_df.stack().reset_index()
        flat_between_df.columns = ["gene_1", "gene_2", "edge"]
        total_between_pairs = flat_between_df.shape[0]

        # count the number of within edges
        num_between_edges = flat_between_df["edge"].sum()
        print(num_between_edges)

        # Get the proportion of 1's looking at within and between genes
        total_pairs = total_within_pairs + total_between_pairs
        prop_within = num_within_edges / total_pairs
        prop_between = num_between_edges / total_pairs
        print("within", prop_within)
        print("between", prop_between)

        # Make output df
        rows.append(
            {
                "Regulon": regulon_name,
                "Lengths": regulon_df.loc[regulon_name, "Lengths"],
                "Genes": geneset,
                "% within edges": prop_within,
                "% between edges": prop_between,
            }
        )
    output_df = pd.DataFrame(rows)
    print(output_df)

    return output_df


pao1_counts_stats = compare_within_btwn_edge(pao1_counts_adj, regulon_df)

pao1_spell_stats = compare_within_btwn_edge(pao1_spell_adj, regulon_df)

print(pao1_counts_stats.shape)
pao1_counts_stats

print(pao1_spell_stats.shape)
pao1_spell_stats.head()

# Save dataframes with threshold info
pao1_counts_stats.to_csv(f"pao1_counts_network_top_{top_percent}.tsv", sep="\t")
pao1_spell_stats.to_csv(f"pao1_spell_network_top_{top_percent}.tsv", sep="\t")

# +
# Format data for plotting
pao1_counts_stats_melt = pd.melt(
    pao1_counts_stats.reset_index(),
    id_vars=["Regulon", "Lengths", "Genes"],
    value_vars=["% within edges", "% between edges"],
)

pao1_spell_stats_melt = pd.melt(
    pao1_spell_stats.reset_index(),
    id_vars=["Regulon", "Lengths", "Genes"],
    value_vars=["% within edges", "% between edges"],
)
# -

print(pao1_counts_stats_melt.shape)
pao1_counts_stats_melt

print(pao1_spell_stats_melt.shape)
pao1_spell_stats_melt

# +
# Make breakdown plot per regulon to see how each contributes
fig1, axes1 = plt.subplots(ncols=1, nrows=2, figsize=(12, 15))

fig1 = sns.barplot(
    data=pao1_counts_stats_melt,
    x="Regulon",
    y="value",
    hue="variable",
    palette="Blues",
    ax=axes1[0],
)
fig1 = sns.barplot(
    data=pao1_spell_stats_melt,
    x="Regulon",
    y="value",
    hue="variable",
    palette="Blues",
    ax=axes1[1],
)
plt.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3
)
axes1[0].set_xticklabels(axes1[0].get_xticklabels(), rotation=40, ha="right")
axes1[1].set_xticklabels(axes1[1].get_xticklabels(), rotation=40, ha="right")
axes1[0].set_title("counts matrix")
axes1[1].set_title("SPELL matrix")
axes1[0].set_ylabel("% edges")
axes1[1].set_ylabel("% edges")
axes1[0].set_xlabel("")
axes1[1].set_xlabel("Regulon")
plt.suptitle(
    f"% within vs between edges (using top {top_percent}%) per regulon", fontsize=14
)
# -

# **Check:**
#
# We noticed that there is a spike in the % of between edges for the *Anr_short_list*. We suspect that these between genes are connected to the other Anr_regulon genes that were removed in this short list. Let's check this.

pao1_spell_stats_test = pao1_spell_stats.set_index("Regulon")

pao1_spell_stats_test.loc[["Anr_regulon", "Anr_short_list"]]

# Get genes from the Anr regulon and the short list
anr_regulon_genes = pao1_spell_stats_test.loc["Anr_regulon", "Genes"]
anr_short_genes = pao1_spell_stats_test.loc["Anr_short_list", "Genes"]

# Get genes removed from the short list
anr_short_missing_genes = list(set(anr_regulon_genes).difference(set(anr_short_genes)))
print(len(anr_short_missing_genes))
assert len(anr_short_missing_genes) == 72 - 5

# +
# Examine between relationships
# Get within edges
within_df = pao1_spell_adj.loc[anr_short_genes, anr_short_genes]
tril_within_df = within_df.where(~np.triu(np.ones(within_df.shape)).astype(np.bool))
flat_within_df = tril_within_df.stack().reset_index()
flat_within_df.columns = ["gene_1", "gene_2", "edge"]
total_within_pairs = flat_within_df.shape[0]

# count the number of within edges
num_within_edges = flat_within_df["edge"].sum()
print(num_within_edges)
# -

pao1_spell_adj.loc[anr_short_genes, anr_short_missing_genes]

# +
# Get between edges
between_df = pao1_spell_adj.loc[anr_short_genes, anr_short_missing_genes]
tril_between_df = between_df.where(~np.triu(np.ones(between_df.shape)).astype(np.bool))
flat_between_df = tril_between_df.stack().reset_index()
flat_between_df.columns = ["gene_1", "gene_2", "edge"]
total_between_pairs = flat_between_df.shape[0]

# count the number of within edges
num_between_edges = flat_between_df["edge"].sum()
print(num_between_edges)
# -

# Get the proportion of 1's looking at within and between genes
total_pairs = total_within_pairs + total_between_pairs
prop_within = num_within_edges / total_pairs
prop_between = num_between_edges / total_pairs
print(prop_within, prop_between)

# +
# Make boxplot for number of edges within vs between genes in gene sets/regulons
fig2, axes2 = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))

fig2 = sns.boxplot(
    data=pao1_counts_stats_melt,
    x="variable",
    y="value",
    palette="Blues",
    notch=True,
    ax=axes2[0],
)

fig2 = sns.boxplot(
    data=pao1_spell_stats_melt,
    x="variable",
    y="value",
    palette="Blues",
    notch=True,
    ax=axes2[1],
)
axes2[0].set_title("counts matrix")
axes2[1].set_title("SPELL matrix")
axes2[0].set_ylabel("% edges")
axes2[1].set_ylabel("% edges")
axes2[0].set_xlabel("")
axes2[1].set_xlabel("")
plt.suptitle(f"% within vs between edges (using top {top_percent}%)", fontsize=14)

# +
# Save plot
# fig1.figure.savefig(f"pao1_within_vs_between_edges_top_{top_percent}_per_regulon.svg", dpi=300)
# fig2.figure.savefig(f"pao1_within_vs_between_edges_top_{top_percent}.svg", dpi=300)
