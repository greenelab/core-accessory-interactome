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
#
# For this experiment, given a set of known gene regulons, that we expect to be clustered together, let's we evaluate which correlation matrix captures these geneset relationships better.
# To determine this we calculate the percentage of within (i.e. edges connecting two genes within the geneset) compared to the percentage of edges between (i.e.edges connecting a gene within the geneset and some other external gene).

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import paths

# ## User params

# +
# Subset of genes
subset_genes = "all"

# Threshold to use to define edges between genes
# Top X% of genes are used
top_percent = 0.01

# Regulon file
regulon_filename = "gene_sets_refs.csv"
# regulon_filename = "regprecise_format.txt"

regulon_delim = ","
# -

# ## Load correlation matrix

pao1_corr_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_raw_mat_test.tsv"
)
pao1_corr_spell_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_log_spell_mat_test.tsv"
)

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

# ### Quick check
# Look at what the distribution of correlation scores looks like to make sure the threshold makes sense

# +
tril_corr_counts_test = pao1_corr_counts.where(
    ~np.triu(np.ones(pao1_corr_counts.shape)).astype(np.bool)
)
tril_corr_spell_test = pao1_corr_spell.where(
    ~np.triu(np.ones(pao1_corr_spell.shape)).astype(np.bool)
)

# Flatten dataframe
flat_corr_counts_test = tril_corr_counts_test.stack().reset_index()
flat_corr_counts_test.columns = ["gene_1", "gene_2", "corr_value"]

flat_corr_spell_test = tril_corr_spell_test.stack().reset_index()
flat_corr_spell_test.columns = ["gene_1", "gene_2", "corr_value"]
# -

sns.displot(flat_corr_counts_test["corr_value"])
sns.displot(flat_corr_spell_test["corr_value"])

# We will create adjacency matrix using threshold defined above.
# The adjacency matrix will determine the strength of the connection between two genes.
# If the concordance is strong enough (i.e. above the threshold), then the genes are connected by an edge.
#
# Since we're examining the connectedness of genes within a regulon we will **not** take the absolute value of the correlation scores. Regulons are groups of operons that share the same transcription factor and therefore we'd expect genes within the same operon to be positively correlated. Operons are a groups of contiguous genes with a shared promoter sequence. For more information the difference between regulon vs operon read [here](https://pediaa.com/what-is-the-difference-between-operon-and-regulon/)

pao1_counts_adj = (pao1_corr_counts > pao1_corr_counts_threshold).astype(float)
pao1_spell_adj = (pao1_corr_spell > pao1_corr_spell_threshold).astype(float)

pao1_counts_adj.head()

pao1_spell_adj.head()

# ## Evaluate relationships captured
#
# For this experiment, given a set of known gene regulons, that we expect to be clustered together, let's we evaluate which correlation matrix captures these geneset relationships better.
# To determine this we calculate the percentage of within (i.e. edges connecting two genes within the geneset) compared to the percentage of edges between (i.e.edges connecting a gene within the geneset and some other external gene).

# Load regulon data
regulon_df = pd.read_csv(regulon_filename, sep=regulon_delim, header=0, index_col=0)

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

        # Processing
        # Since there are some gene ids that are from PA14
        # We will take the intersection
        geneset_processed = set(geneset).intersection(adj_df.index)

        # WITHIN gene calculations
        # Get within edges
        # Since this is a symmetric matrix we only will consider the
        # bottom triangle (excluding the diagonal) in order to count
        # the number of unique connections
        within_df = adj_df.loc[geneset_processed, geneset_processed]
        tril_within_df = within_df.where(
            ~np.triu(np.ones(within_df.shape)).astype(np.bool)
        )

        flat_within_df = tril_within_df.stack().reset_index()
        flat_within_df.columns = ["gene_1", "gene_2", "edge"]
        total_within_pairs = flat_within_df.shape[0]

        # Get proportion of within edges
        num_within_edges = flat_within_df["edge"].sum()
        prop_within = num_within_edges / total_within_pairs

        # BETWEEN calculations
        # Get between edges
        # This matrix is a rectangle of genes in the regulon x genes not in the regulon
        # so we can look at connectsion across the entire matrix as opposed to taking
        # the lower triangle
        not_geneset = set(adj_df.index).difference(geneset_processed)
        between_df = adj_df.loc[geneset_processed, not_geneset]

        flat_between_df = between_df.stack().reset_index()
        flat_between_df.columns = ["gene_1", "gene_2", "edge"]
        total_between_pairs = flat_between_df.shape[0]

        # count the number of within edges
        num_between_edges = flat_between_df["edge"].sum()
        prop_between = num_between_edges / total_between_pairs

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
axes1[0].set_title(f"counts matrix (threshold={pao1_corr_counts_threshold: .2f})")
axes1[1].set_title(f"SPELL matrix (threshold={pao1_corr_spell_threshold: .2f})")
axes1[0].set_ylabel("% edges")
axes1[1].set_ylabel("% edges")
axes1[0].set_xlabel("")
axes1[1].set_xlabel("Regulon")
plt.suptitle(
    f"% within vs between edges (using top {top_percent*100}%) per regulon", fontsize=14
)

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
axes2[0].set_title(f"counts matrix (threshold={pao1_corr_counts_threshold: .2f})")
axes2[1].set_title(f"SPELL matrix (threshold={pao1_corr_spell_threshold: .2f})")
axes2[0].set_ylabel("% edges")
axes2[1].set_ylabel("% edges")
axes2[0].set_xlabel("")
axes2[1].set_xlabel("")
plt.suptitle(f"% within vs between edges (using top {top_percent*100}%)", fontsize=14)
# -

# Look at the regulons that differed the most to make sure there isn't some sort of bias -
# like they all have very few genes
# Doesn't look to be the case
regulon_df.loc[
    [
        "Zur_regulon",
        "Anr_short_list",
        "PhoB_short_list",
        "AlgU_short_list",
        "LasR_short_list",
        "RhlR_short_list",
        "GbdR_regulon",
        "ErdR_regulon",
        "SoxR_regulon",
        "PhhR_regulon",
    ]
]

# Save plot
fig1.figure.savefig(
    f"pao1_within_vs_between_edges_top_{top_percent}_per_regulon.svg", dpi=300
)
fig2.figure.savefig(f"pao1_within_vs_between_edges_top_{top_percent}.svg", dpi=300)


# ### Check
#
# As a control, we will create pseudo-regulons by creating regulons with a random set of size matched genes and comparing the within versus between connections.
#
# The way we calculate the proportions of within and between connections is to take the number of within connections divided by the number of all possible within pairs.
# Similarly for the between connections we take the number of between connections divided by the number of all possible between pairs. For this calculation, the denominators are separate because the number of non-regulon genes is much larger.
# Our _hypothesis_ is that the between connections are mostly due to noise and randomness in the data. If we divided the number of within connections by the total number of all pairs (within pairs + between pairs) then any signal from the within connections will get drowned out by the noise from the between connections because the size of the regulons are so much smaller.
#
# If the between connections are mostly noise, then we'd expect using a shuffled set of regulons, our within and between connections to be roughly equal in distribution. This is what we find.

# Make shuffled regulon dataset
def make_shuffled_regulon(regulon_df, gene_ids):
    """
    This function inputs real regulon data and creates
    a fake regulon dataset by randomly sampling groups of
    size matched regulons
    """
    rows = []
    for regulon_name in regulon_df.index:
        len_regulon = regulon_df.loc[regulon_name, "Lengths"]

        # Select random set of size-matched genes
        random_geneset = random.sample(gene_ids, len_regulon)

        # Make output df
        rows.append(
            {
                "Regulon": regulon_name,
                "Lengths": regulon_df.loc[regulon_name, "Lengths"],
                "Genes": random_geneset,
            }
        )
    output_df = pd.DataFrame(rows)

    return output_df


shuffled_regulon_df = make_shuffled_regulon(regulon_df, list(pao1_corr_counts.index))

shuffled_regulon_df.head()

pao1_counts_shuffled_stats = compare_within_btwn_edge(
    pao1_counts_adj, shuffled_regulon_df
)

pao1_spell_shuffled_stats = compare_within_btwn_edge(
    pao1_spell_adj, shuffled_regulon_df
)

pao1_counts_shuffled_stats.head()

# +
# Format data for plotting
pao1_counts_shuffled_stats_melt = pd.melt(
    pao1_counts_shuffled_stats.reset_index(),
    id_vars=["Regulon", "Lengths", "Genes"],
    value_vars=["% within edges", "% between edges"],
)

pao1_spell_shuffled_stats_melt = pd.melt(
    pao1_spell_shuffled_stats.reset_index(),
    id_vars=["Regulon", "Lengths", "Genes"],
    value_vars=["% within edges", "% between edges"],
)
# -

pao1_counts_shuffled_stats_melt.head()

# +
# Make boxplot for number of edges within vs between genes in gene sets/regulons
fig3, axes3 = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))

fig3 = sns.boxplot(
    data=pao1_counts_shuffled_stats_melt,
    x="variable",
    y="value",
    palette="Blues",
    notch=True,
    ax=axes3[0],
)

fig3 = sns.boxplot(
    data=pao1_spell_shuffled_stats_melt,
    x="variable",
    y="value",
    palette="Blues",
    notch=True,
    ax=axes3[1],
)
axes3[0].set_title(f"counts matrix (threshold={pao1_corr_counts_threshold: .2f})")
axes3[1].set_title(f"SPELL matrix (threshold={pao1_corr_spell_threshold: .2f})")
axes3[0].set_ylabel("% edges")
axes3[1].set_ylabel("% edges")
axes3[0].set_xlabel("")
axes3[1].set_xlabel("")
plt.suptitle(
    f"% SHUFFLED within vs between edges (using top {top_percent*100}%)", fontsize=14
)
