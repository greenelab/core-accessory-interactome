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

# +
# Save
# pao1_adj_filename = os.path.join(paths.LOCAL_DATA_DIR, f"pao1_adj_threshold_{corr_threshold}.tsv")
# pao1_adj.to_csv(pao1_adj_filename, sep="\t")
# -

# ## Module detection
# Detect modules. Get membership of genes that are closely related based on adjacency matrix (Function: `TOMsimilarity`)
#
# * First we need to calculate the topological overlap measure (TOM). The topological overlap of two nodes reflects their similarity in terms of the commonality of the nodes they connect to. What is the input and output???
#
# * We can group genes into modules based on their shared connections using clustering. What is the input and output of clustering???

# + language="R"
# library("WGCNA")

# + magic_args="-i pao1_adj -o TOMdis" language="R"
#
# pao1_adj_mat <- as.matrix(pao1_adj)
# print(is.numeric(pao1_adj_mat))
#
# # Similarity based on adjacency
# TOMsim <- TOMsimilarity(pao1_adj_mat)
# TOMdis <- 1-TOMsim
#
# mode(TOMdis) <- "integer"
#
# # clustering
# modules <- hclust(TOMdis)

# +
# TO DO
# cutTreeDynamic to get modules from htree?

# Format output
# Get community membership for a single threshold
# Format and save output to have columns: gene_id | group_id

# Organize and clean
# Organize script to be able to run for multiple different thresholds as a script and save membership
# Add comments on what each step is doing
