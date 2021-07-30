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
# This notebook performs community detection approaches to identify network modules. Community detection considers genes to be in a network/graph where genes are connected to other genes based on similarity between expression pattern across samples (i.e. correlation score between gene A and B). Community detection will define modules as a group of genes that are densely connected to each other but sparsely connected to other genes in the network (within vs between edges). Here each gene still belongs to a single module
#
# The output of this notebook are files that have the following columns:
# gene id | module id
#
# Note: All methods here are using undirected weighted networks. All methods take edge weights as input.

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
method = "infomap"

# Params for different methods to adjust
# length of random walk to perform for walktrap
# Short random walks tend to stay in the same community
nsteps = 10

# Number of trials to attempt to partition the network for infomap
ntrials = 20

# TO DO
# Think about running these methods across multiple seeds and taking a consensus

# +
# Load correlation matrix
pao1_pearson_mat_filename = paths.PAO1_CORR_LOG_SPELL
pa14_pearson_mat_filename = paths.PA14_CORR_LOG_SPELL

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

print(pao1_corr_graph.shape)
pao1_corr_graph.head()

# Make into a graph object
pao1_G = ig.Graph.TupleList(pao1_corr_graph.values, weights=True, directed=False)
pa14_G = ig.Graph.TupleList(pa14_corr_graph.values, weights=True, directed=False)

# make sure vertex/edge properties exist
print(pao1_G.es["weight"][:5])

# ## Community detection

# ### Fast-greedy
# This algorithm starts from a completely unclustered set of nodes and iteratively adds communities such that the modularity (score maximizing within edges and minimizing between edges) is maximized until no additional improvement can be made.
#
# Note: Looks like fast-greedy requires a simple graph (i.e. no multiple edges per node), so we use [simplify](https://igraph.org/python/doc/api/igraph._igraph.GraphBase.html#simplify) to combine edges
#
# Returns VertextDendrogram

# %%time
if method == "fastgreedy":
    # Simplify graph to remove multiple edges and loops
    if not pao1_G.is_simple():
        pao1_G.simplify()
    if not pa14_G.is_simple():
        pa14_G.simplify()

    assert pao1_G.is_simple()
    assert pa14_G.is_simple()

    # Detection method
    pao1_partition = pao1_G.community_fastgreedy(weights=pao1_G.es["weight"])
    pa14_partition = pa14_G.community_fastgreedy(weights=pa14_G.es["weight"])

# +
# pao1_G.vs.attribute_names()
# -

# ### Walktrap
# This algorithm performs random walks using a specified step size. Where densely connected areas occur, the random walk becomes “trapped” in local regions that then define communities
#
# Returns VertextDendrogram

# %%time
if method == "walktrap":
    pao1_partition = pao1_G.community_walktrap(
        weights=pao1_G.es["weight"], steps=nsteps
    )
    pa14_partition = pa14_G.community_walktrap(
        weights=pa14_G.es["weight"], steps=nsteps
    )

# ### Multilevel
# This algorithm is similar to fastgreedy, but it merges communities to optimize modularity based upon only the neighboring communities as opposed to all communities. The algorithm terminates when only a single node is left, or when the improvement in modularity cannot result from the simple merge of two neighboring communities. (Louvain clustering)
#
# Returns VertexClustering

# %%time
if method == "louvain":
    pao1_partition = pao1_G.community_multilevel(
        weights=pao1_G.es["weight"], return_levels=False
    )
    pa14_partition = pa14_G.community_multilevel(
        weights=pa14_G.es["weight"], return_levels=False
    )

# ### Infomap
# This algorithm uses the probability flow of information in random walks, which occurs more readily in groups of heavily connected nodes. Thus, information about network structure can be compressed in maps of modules (nodes where information travels quickly)
#
# Returns VertexClustering

# %%time
if method == "infomap":
    pao1_partition = pao1_G.community_infomap(
        edge_weights=pao1_G.es["weight"], trials=ntrials
    )
    pa14_partition = pa14_G.community_infomap(
        edge_weights=pa14_G.es["weight"], trials=ntrials
    )


# ## Get membership

# get dataframe mapping Pa genes to communities
def graph_partition_to_df(G, partition, method):
    if method in ["louvain", "infomap"]:
        clusters = []
        for label, vl in enumerate(partition):
            clusters += [(G.vs["name"][v], label, G.degree(v)) for v in vl]

        membership_df = pd.DataFrame(clusters, columns=["gene", "module id", "degree"])
        membership_df = membership_df.set_index("gene")

        return membership_df


# +
# pao1_partition.es.attribute_names()
# -

pao1_membership_df = graph_partition_to_df(pao1_G, pao1_partition, method)
print(len(pao1_membership_df["module id"].unique()))
pao1_membership_df.sort_values(by="degree", ascending=False).head()

pa14_membership_df = graph_partition_to_df(pa14_G, pa14_partition, method)
print(len(pa14_membership_df["module id"].unique()))
pa14_membership_df.sort_values(by="degree", ascending=False).head()

# Save membership dataframe
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{method}.tsv"
)
pao1_membership_df.to_csv(pao1_membership_filename, sep="\t")
pa14_membership_df.to_csv(pa14_membership_filename, sep="\t")
