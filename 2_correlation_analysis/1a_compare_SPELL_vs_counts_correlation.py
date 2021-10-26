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
import os
import pandas as pd
import numpy as np
from scripts import paths

# ## User params

# +
# Which correlation matrix to use. Choices = ["raw", "spell"]
which_corr = "raw"

# Threshold to use to define edges between genes
# Top X% of genes are used
top_percent = 0.1
# -

# ## Load correlation matrix

if which_corr == "raw":
    pao1_corr_filename = paths.PAO1_CORR_RAW
    pa14_corr_filename = paths.PA14_CORR_RAW
elif which_corr == "spell":
    pao1_corr_filename = paths.PAO1_CORR_LOG_SPELL
    pa14_corr_filename = paths.PA14_CORR_LOG_SPELL

# Load correlation data
pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pa14_corr = pd.read_csv(pa14_corr_filename, sep="\t", index_col=0, header=0)


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


pao1_corr_threshold = get_corr_threshold(pao1_corr, top_percent)

pa14_corr_threshold = get_corr_threshold(pa14_corr, top_percent)

# Create adjacency matrix using threshold defined above
# The adjacency matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by an edge
# TO DO:abs?????
pao1_adj = (pao1_corr > pao1_corr_threshold).astype(float)
pa14_adj = (pa14_corr > pa14_corr_threshold).astype(float)

pao1_adj.head()

# ## Examine relationships
#
# Given a set of known gene regulons, that we expect to be clustered together, let's calculate if that is the case.
#
# Here we calculate the percentage of within (i.e. edges connecting genes from the s vs between edges there are between genes from a given regulon.

regulon_filename = "gene_sets_refs.csv"

regulon_df = pd.read_csv(regulon_filename, header=0, index_col=0)

regulon_df["Genes"] = regulon_df["Genes"].str.split(";").apply(list)

regulon_df


def calc_within_edge(adj_df, regulon_df):
    # Loop through each regulon
    for regulon_name in regulon_df.index:
        geneset = regulon_df.loc[regulon_name, "Genes"]
        print(geneset)
        print(len(geneset))

        # Since there are some gene ids that are from PA14
        # We will take the intersection
        geneset_processed = set(geneset).intersection(adj_df.index)

        # Get within edges
        within_genes = adj_df.loc[geneset_processed, geneset_processed]
        print(within_genes)

        # Get between edges
        not_geneset = set(adj_df.index).difference(geneset_processed)
        between_genes = adj_df.loc[geneset_processed, not_geneset]

        print(between_genes)

        # Get the proportion of 1's looking at within and between genes

        # make df

        break

        # return


calc_within_edge(pao1_adj, regulon_df)

# +
# Make boxplot for number of edges within vs between genes in gene sets/regulons

# +
# Repeat this for different thresholds and save each one (spell, counts) (0.01, 0.05, 0.1, 0.2)

# +
# Try this for single threshold for one dataset for other gene sets
# -

"""# Save membership dataframe
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{cluster_method}_{gene_subset}_{processed}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{cluster_method}_{gene_subset}_{processed}.tsv"
)
pao1_membership_df.to_csv(pao1_membership_filename, sep="\t")
pa14_membership_df.to_csv(pa14_membership_filename, sep="\t")"""
