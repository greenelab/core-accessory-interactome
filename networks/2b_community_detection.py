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

# # Network community detection
#
# This notebook performs community detection approaches to identify network modules.

# +
import os
import random

import numpy as np
import pandas as pd
import igraph as ig
from core_acc_modules import paths

# +
# User params

# Choices = ["fastgreedy", "walktrap", "louvain", "infomap"]
method = "fastgreedy"

# TO DO: params for different methods to adjust
# steps for walktrap
# trails for infomap

# +
# Load correlation matrix --> which correlation matrix to use?
pao1_pearson_mat_filename = os.path.join(paths.LOCAL_DATA_DIR, "pao1_pearson_mat.tsv")
pa14_pearson_mat_filename = os.path.join(paths.LOCAL_DATA_DIR, "pa14_pearson_mat.tsv")

# Take abs of correlation scores
# In this case we care about the strength and not the direction
pao1_corr = pd.read_csv(
    pao1_pearson_mat_filename, sep="\t", index_col=0, header=0
).abs()
pa14_corr = pd.read_csv(
    pa14_pearson_mat_filename, sep="\t", index_col=0, header=0
).abs()
# -

pao1_corr.head()

# +
# Format correlation matrix into graph (i.e. dataframe with edge weight per pair of genes)
# The dataframe should have columns: from, to, weight
pao1_corr_graph = pao1_corr.stack().reset_index()
pao1_corr_graph.columns = ["from", "to", "weight"]

pa14_corr_graph = pa14_corr.stack().reset_index()
pa14_corr_graph.columns = ["from", "to", "weight"]
# -

# Drop duplicate rows since correlation matrix is symmetric
pao1_corr_graph = pao1_corr_graph.drop_duplicates()
pa14_corr_graph = pa14_corr_graph.drop_duplicates()

# Drop gene loops
# Note 'query' not working for some reason
pao1_corr_graph = pao1_corr_graph[pao1_corr_graph["from"] != pao1_corr_graph["to"]]
pa14_corr_graph = pa14_corr_graph[pa14_corr_graph["from"] != pa14_corr_graph["to"]]

pao1_corr_graph.head()

# Make into a graph object
pao1_G = ig.Graph.TupleList(pao1_corr_graph.values, weights=True, directed=False)
pa14_G = ig.Graph.TupleList(pa14_corr_graph.values, weights=True, directed=False)

# make sure vertex/edge properties exist
print(pao1_G.es["weight"][:5])

# +
# TO DO: Add label for core, accessory gene
# -

# ## Community detection

# ### Fast-greedy
# This algorithm starts from a completely unclustered set of nodes and iteratively adds communities such that the modularity (score maximizing within edges and minimizing between edges) is maximized until no additional improvement can be made.
#
# **What is this simplification step doing?**
# This is removing multiple edges and loops -- how???

if method == "fastgreedy":
    pao1_partition = pao1_G.simplify().community_fastgreedy(weights=pao1_G.es["weight"])
    pa14_partition = pao1_G.simplify().community_fastgreedy(weights=pa14_G.es["weight"])

# +
# Error at fast_community.c:553: fast-greedy community finding works only on graphs without multiple edges, Invalid value
# -

# ### Walktrap
# This algorithm performs random walks using a specified step size. Where densely connected areas occur, the random walk becomes “trapped” in local regions that then define communities
#

if method == "walktrap":
    pao1_partition = pao1_G.community_walktrap(weights=pao1_G.es["weight"])
    pa14_partition = pa14_G.community_walktrap(weights=pao1_G.es["weight"])

# ### Multilevel
# This algorithm is similar to fastgreedy, but it merges communities to optimize modularity based upon only the neighboring communities as opposed to all communities. The algorithm terminates when only a single node is left, or when the improvement in modularity cannot result from the simple merge of two neighboring communities. (Louvain clustering)

if method == "louvain":
    pao1_partition = pao1_G.community_multilevel(
        weights=pao1_G.es["weight"], return_levels=False
    )
    pa14_partition = pa14_G.community_multilevel(
        weights=pao1_G.es["weight"], return_levels=False
    )

# ### Infomap
# This algorithm uses the probability flow of information in random walks, which occurs more readily in groups of heavily connected nodes. Thus, information about network structure can be compressed in maps of modules (nodes where information travels quickly)
#

if method == "infomap":
    pao1_partition = pao1_G.community_infomap(edge_weights=pao1_G.es["weight"])
    pa14_partition = pa14_G.community_infomap(edge_weights=pao1_G.es["weight"])


# ## Get membership

# get dataframe mapping Pa genes to communities
def graph_partition_to_df(G, partition):
    clusters = []
    for label, vl in enumerate(partition):
        clusters += [(G.vs["name"][v], label, G.degree(v)) for v in vl]
    return pd.DataFrame(clusters, columns=["gene", "module id", "degree"])


pao1_membership_df = graph_partition_to_df(pao1_G, pao1_partition)
print(len(pao1_membership_df["module id"].unique()))
pao1_membership_df.sort_values(by="degree", ascending=False).head()

pa14_membership_df = graph_partition_to_df(pa14_G, pa14_partition)
print(len(pa14_membership_df["module id"].unique()))
pa14_membership_df.sort_values(by="degree", ascending=False).head()

# Save
# Save membership dataframe
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{method}.tsv"
)
pao1_membership_df.to_csv(pao1_membership_filename, sep="\t")
pa14_membership_df.to_csv(pa14_membership_filename, sep="\t")
