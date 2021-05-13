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

# # Get network communities
#
# This notebook gets network communities for the compendia (PAO1 and PA14) using different thresholds.
#
# The output of this notebook are files for each threshold. These files have the following columns:
# gene id | module id

# +
# %load_ext autoreload
# %autoreload 2
# %load_ext rpy2.ipython
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from core_acc_modules import paths
from rpy2.robjects import pandas2ri

pandas2ri.activate()
# -

# ## Set user parameters
#
# For now we will vary the correlation threshold (`corr_threshold`) but keep the other parameters consistent
#
# We will run this notebook for each threshold parameter

# +
# Params
corr_threshold = 0.5

# Params for hclust
# clustering_method is the distance metric used to determine if clusters should be merged
# https://en.wikipedia.org/wiki/Hierarchical_clustering
clustering_method = "average"

# Params for cutreeDynamic
# minimum cluster size
min_cluster_size = 30

# The higher the value (or if TRUE), the more and smaller clusters will be produced
deep_split = 2


# Output files
pao1_membership_filename = f"pao1_membership_{corr_threshold}.tsv"
pa14_membership_filename = f"pa14_membership_{corr_threshold}.tsv"
# -

# Load expression data
pao1_compendium_filename = paths.PAO1_COMPENDIUM
pa14_compendium_filename = paths.PA14_COMPENDIUM

pao1_compendium = pd.read_csv(pao1_compendium_filename, sep="\t", header=0, index_col=0)
pa14_compendium = pd.read_csv(pa14_compendium_filename, sep="\t", header=0, index_col=0)

print(pao1_compendium.shape)
pao1_compendium.head()

print(pa14_compendium.shape)
pa14_compendium.head()

# ## Make adjacency matrix

# Get perason correlation
# This correlation matrix will represent the concordance
# between two gene expression profiles
pao1_corr = pao1_compendium.corr()
pa14_corr = pa14_compendium.corr()

pao1_corr.head()

# Create adjacency matrix using threshold defined above
# The adjacency matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by an edge
pao1_adj = (pao1_corr.abs() >= corr_threshold).astype(float)
pa14_adj = (pa14_corr.abs() >= corr_threshold).astype(float)

pao1_adj.head()

# ## Plot
#
# Plot clustering of adjacency matrix

# +
# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_adj, cmap="viridis")
h1.fig.suptitle(f"Adjacency of PAO1 genes using threshold={corr_threshold}")

# Save
pao1_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_adj_{corr_threshold}_clustermap.png"
)
h1.savefig(pao1_clustermap_filename, dpi=300)

# +
# Plot heatmap
plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_adj, cmap="viridis")
h2.fig.suptitle(f"Adjacency of PA14 genes using threshold={corr_threshold}")

# Save
pa14_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_adj_{corr_threshold}_clustermap.png"
)
h2.savefig(pa14_clustermap_filename, dpi=300)
# -

# ## Module detection
# To detect modules, we want to look for genes that are closely related based on the adjacency matrix (i.e. genes that have similar connections)
#
# First we need to calculate the topological overlap measure (TOM) using [TOMsimilarity](https://rdrr.io/cran/WGCNA/man/TOMsimilarity.html). The topological overlap of two nodes reflects their similarity in terms of the commonality of the nodes they connect to.
#
# * input: adjacency matrix (square symmetric matrix with 0 and 1 entries)
# * output: matrix holding the topological overlap. For an unweighted network, topological overlap = 1 if node _i_ and _j_ are linked and the neighbors of the node _i_ is connected to all of the neighbors of node _j_. Topological overlap = 0 if node _i_ and _j_ are unlinked and the two nodes have no common neighbors. The TOM is a matrix of continuous values ranging from 0 to 1.
#
# Next, we need to cluster based on this information. For now, we will use heirarchal clustering, [hclust](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust) to identify communities. A community is a group of genes with shared connections. This function performs a hierarchical cluster analysis using a set of dissimilarities for the n objects being clustered. Initially, each object is assigned to its own cluster and then the algorithm proceeds iteratively, at each stage joining the two most similar clusters, continuing until there is just a single cluster. At each stage distances between clusters are recomputed by the Lance–Williams dissimilarity update formula according to the particular clustering method being used.
#
#
# * input: distance of the dissimilarity matrix (i.e. distance of 1- TOM). The dissimilarity matrix is a matrix that is max = 1 if two genes are very dissimilar based on their connections. This matrix contains continuous values. Taking the distance of the dissimilarity will compare the genes based on their dissimilarity scores and return a distance between them. This distance matrix is a continuous value. This matrix is a lower triangle matrix.
# * output: hclust object which describes the [tree](https://en.wikipedia.org/wiki/Dendrogram) produced by the clustering process.
#
# Finally we need to use [cutreeDynamic](https://rdrr.io/cran/dynamicTreeCut/man/cutreeDynamic.html) to identify modules based on hclust output (tree). Cutting the tree at a given height will give a partitioning clustering at a selected precision.
#
# * input: hierarchial clusterig dendogram returned from hclust
# * output: A vector of numerical labels giving assignment of objects to modules. Unassigned objects are labeled 0, the largest module has label 1, next largest 2 etc.

# + language="R"
# library("WGCNA")

# + magic_args="-i pao1_adj -i clustering_method -i deep_split -i min_cluster_size -o pao1_modules -o TOM_pao1" language="R"
#
# pao1_adj_mat <- as.matrix(pao1_adj)
# print(is.numeric(pao1_adj_mat))
#
# # Similarity based on adjacency
# TOM_pao1 <- TOMsimilarity(pao1_adj_mat)
# dissTOM <- 1-TOM_pao1
#
# # Clustering
# geneTree <- hclust(as.dist(dissTOM), method=clustering_method)
#
# # Module identification using dynamic tree cut:
# pao1_modules <- cutreeDynamic(
#     dendro = geneTree,
#     distM = dissTOM,
#     deepSplit = deep_split,
#     minClusterSize = min_cluster_size
# )
#
#
# table(pao1_modules)

# + magic_args="-i pa14_adj -i clustering_method -i deep_split -i min_cluster_size -o pa14_modules -o TOM_pa14" language="R"
#
# pa14_adj_mat <- as.matrix(pa14_adj)
# print(is.numeric(pa14_adj_mat))
#
# # Similarity based on adjacency
# TOM_pa14 <- TOMsimilarity(pa14_adj_mat)
# dissTOM <- 1-TOM_pa14
#
# print("TOM")
# print(TOM_pa14)
#
# # Clustering
# geneTree <- hclust(as.dist(dissTOM), method=clustering_method)
#
# # Module identification using dynamic tree cut:
# pa14_modules <- cutreeDynamic(
#     dendro = geneTree,
#     distM = dissTOM,
#     deepSplit = deep_split,
#     minClusterSize = min_cluster_size
# )
#
#
# table(pa14_modules)
# -

# ## Plot TOM

# +
# Plot heatmap
plt.figure(figsize=(20, 20))
h3 = sns.clustermap(TOM_pao1, cmap="viridis")
h3.fig.suptitle(f"Similarity (TOM) of PAO1 genes using threshold={corr_threshold}")

# Save
pao1_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_TOM_{corr_threshold}_clustermap.png"
)
h3.savefig(pao1_clustermap_filename, dpi=300)

# +
# Plot heatmap
plt.figure(figsize=(20, 20))
h4 = sns.clustermap(TOM_pa14, cmap="viridis")
h4.fig.suptitle(f"Similarity (TOM) of PA14 genes using threshold={corr_threshold}")

# Save
pa14_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_TOM_{corr_threshold}_clustermap.png"
)
h4.savefig(pa14_clustermap_filename, dpi=300)
# -

# ## Get membership

# +
# Get module membership for a single threshold
# Format and save output to have columns: gene_id | group_id
pao1_membership_df = pd.DataFrame(
    data={"module id": pao1_modules}, index=pao1_adj.index
)

pao1_membership_df["module id"].value_counts()

# +
# Get module membership for a single threshold
# Format and save output to have columns: gene_id | group_id
pa14_membership_df = pd.DataFrame(
    data={"module id": pa14_modules}, index=pa14_adj.index
)

pa14_membership_df["module id"].value_counts()

# +
# Save membership dataframe
# pao1_membership_df.to_csv(pao1_membership_filename, sep="\t")
# pa14_membership_df.to_csv(pa14_membership_filename, sep="\t")
