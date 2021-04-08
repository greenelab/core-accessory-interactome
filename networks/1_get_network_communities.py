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
#     display_name: Python [conda env:core_acc_env] *
#     language: python
#     name: conda-env-core_acc_env-py
# ---

# # Get network communities
#
# This notebook gets network communities for the compendia (PAO1 and PA14) using different thresholds.
#
# The output of this notebook are files for each threshold. These files have the following columns:
# gene_id | group_id

# +
# %load_ext autoreload
# %autoreload 2
# %load_ext rpy2.ipython
import os
import pandas as pd
from core_acc_modules import paths
from rpy2.robjects import pandas2ri

pandas2ri.activate()
# -

# Params
corr_threshold = 0.9

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
#
# * Start with similarity matrix, S = [s_ij] = |cor(i,j)| \in [0,1] which represents the concordance between gene expression profiles for gene i and j.
# * The adjacency matrix, A = [a_ij]  encodes the strength of the connection between gene i and j.
# * For 'hard thresholding' a_ij = 1 if s_ij >= threshold, else 0.
# * For 'soft thresholding' a_ij = |s_ij|^\beta . --
#     * Beta is chosen based on 2 criteria.
#     * 1) To get a network follows a scale free topology
#     * Scale-free topology is defined by a network where the probability that a node is connected with k other node (the degree distribution p(k) of a network) decays as a power law p(k) ~ k^(-x)
#     * This is what the “Weighted” in the name refers to -- using a weighted adjacency matrix, where the co-expression similarity is raised to a power
#     * 2) Maintain a high number of connections
# * Here we’re using a hard threshold cutoff
# * We can select different correlation function, default is Pearson.

# Get perason correlation
pao1_corr = pao1_compendium.corr()
pa14_corr = pa14_compendium.corr()

pao1_corr.head()

# Create adjacency matrix using threshold defined above
pao1_adj = (pao1_corr >= corr_threshold).astype(int)
pa14_adj = (pa14_corr >= corr_threshold).astype(int)

pao1_adj.head()

# ## Module detection
# Detect modules. Get membership of genes that are closely related based on adjacency matrix (Function: `TOMsimilarity`)
# * First we need to calculate the topological overlap measure (TOM)
# * The topological overlap of two nodes reflects their similarity in terms of the commonality of the nodes they connect to
# * We can group genes into modules based on their shared connections using clustering -- we’ll probably start with hierarchical clustering

# + language="R"
# library("WGCNA")

# + magic_args="-i pao1_adj -o tom_similarity" language="R"
#
# print(head(pao1_adj))
# tom_similarity = TOMsimilarity(pao1_adj)

# +
# TO DO
# # 1-TOMsim?
# # hclust(1-TOMsim)?
# cutTreeDynamic to get modules from htree?
# Get community membership for a single threshold
# Format and save output to have columns: gene_id | group_id
# Organize script to be able to run for multiple different thresholds as a script and save membership
