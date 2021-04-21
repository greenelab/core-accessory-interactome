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
# gene id | module id

# %load_ext autoreload
# %autoreload 2
import os
import pandas as pd
from sklearn.cluster import DBSCAN, AgglomerativeClustering, AffinityPropagation
from core_acc_modules import paths

# ## Set user parameters
#
# For now we will vary the correlation threshold (`corr_threshold`) but keep the other parameters consistent
#
# We will run this notebook for each threshold parameter

# +
# Params
corr_threshold = 0.9

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

# ## Get similarity between genes
#
# To determine if genes are similar, we will calculate the correlation between genes and apply our threshold cutoff. When we apply our threshold, any scores that are below the threshold will be set to 0 (i.e. using a threshold of 0.5 means that gene pairs that have correlation scores of 0.52 and 0.78 will be left as is but a gene pair with a correlation score of 0.48 will be set to 0).

# Get perason correlation
# This correlation matrix will represent the concordance
# between two gene expression profiles
pao1_corr = pao1_compendium.corr()
pa14_corr = pa14_compendium.corr()

pao1_corr.head()

# +
# Create a similarity matrix usingn the threshold defined above
# The similarity matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by by the correlation score, otherwise the value is set to 0
pao1_corr[pao1_corr.abs() < corr_threshold] = 0.0
pa14_corr[pa14_corr.abs() < corr_threshold] = 0.0

pao1_corr.head()
# -

pa14_corr.head()

# ## Module detection
# To detect modules, will use a clustering algorithm
#
# **About the methods**
#
# * [Hierachal clustering](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html#sklearn.cluster.AgglomerativeClustering):
# * [DBSCAN](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html#sklearn.cluster.DBSCAN):
# * [Affinity propogation](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AffinityPropagation.html#sklearn.cluster.AffinityPropagation):
#
# **ADD RANDOM SEEDS, PLAY WITH PARAMS**
#
# What method will we use and why??

# +
# Clustering using DBSCAN
# pao1_clustering = DBSCAN().fit(pao1_corr)
# pa14_clustering = DBSCAN().fit(pa14_corr)

# Clustering using hierarchal clustering
pao1_clustering = AgglomerativeClustering(linkage="average").fit(pao1_corr)
pa14_clustering = AgglomerativeClustering(linkage="average").fit(pa14_corr)

# Clustering using affinity propogation
# pao1_clustering = AffinityPropagation().fit(pao1_corr)
# pa14_clustering = AffinityPropagation().fit(pa14_corr)

# +
# Get module membership for a single threshold
# Format and save output to have columns: gene_id | group_id
pao1_membership_df = pd.DataFrame(
    data={"module id": pao1_clustering.labels_}, index=pao1_corr.index
)

pao1_membership_df["module id"].value_counts()

# +
# Get module membership for a single threshold
# Format and save output to have columns: gene_id | group_id
pa14_membership_df = pd.DataFrame(
    data={"module id": pa14_clustering.labels_}, index=pa14_corr.index
)

pa14_membership_df["module id"].value_counts()

# +
# Save membership dataframe
# pao1_membership_df.to_csv(pao1_membership_filename, sep="\t")
# pa14_membership_df.to_csv(pa14_membership_filename, sep="\t")
