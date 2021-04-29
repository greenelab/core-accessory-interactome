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
corr_threshold = 0.5

# Correlation matrix files
pao1_corr_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_corr_{corr_threshold}.tsv"
)
pa14_corr_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_corr_{corr_threshold}.tsv"
)

# +
## TO DO
# Create variable to determine which clustering method to use
# -

# Load correlation data
pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pa14_corr = pd.read_csv(pa14_corr_filename, sep="\t", index_col=0, header=0)

# ## Module detection
# To detect modules, we will use a clustering algorithm

# ### DBSCAN
# [DBSCAN](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html#sklearn.cluster.DBSCAN):  Density-Based Spatial Clustering of Applications with Noise views clusters as areas of high density separated by areas of low density. The central component to the DBSCAN is the concept of _core samples_, which are samples that are in areas of high density. A cluster is therefore a set of _core samples_ that are close to each other (measured by some distance measure) and a set of non-core samples that are close to a core sample (but are not themselves core samples).
#
# A cluster is a set of core samples that can be built by recursively taking a core sample, finding all of its neighbors that are core samples, finding all of their neighbors that are core samples, and so on. A cluster also has a set of non-core samples, which are samples that are neighbors of a core sample in the cluster but are not themselves core samples. Intuitively, these samples are on the fringes of a cluster.
#
# * We define a core sample as being a sample in the dataset such that there exist `min_samples` other samples within a distance of `eps`, which are defined as neighbors of the core sample.
# * Using all default values
#

# +
# Clustering using DBSCAN
# pao1_clustering = DBSCAN().fit(pao1_corr)
# pa14_clustering = DBSCAN().fit(pa14_corr)
# -

# ### Hierarchical clustering
# [Hierarchical clustering](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html#sklearn.cluster.AgglomerativeClustering): Initially, each object is assigned to its own cluster and then the algorithm proceeds iteratively, at each stage joining the two most similar clusters (i.e. linkage distance is minimized), continuing until there is just a single cluster.
#
# * n_cluster: The number of clusters to find.
# * linkage: Criterion used to determine distance between observations. 'average'=average distance of each observation in the two sets.
#
# * Note: It looks like this method tends to produce 1 very large cluster. To break this up we will iteratively apply hierarchal clustering on the largest cluster.

# Clustering using hierarchal clustering
pao1_clustering = AgglomerativeClustering(
    n_clusters=None, distance_threshold=0.5, linkage="average"
).fit(pao1_corr)
pa14_clustering = AgglomerativeClustering(
    n_clusters=None, distance_threshold=0.5, linkage="average"
).fit(pa14_corr)

# ### Affinity propogation
#
# [Affinity propogation](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AffinityPropagation.html#sklearn.cluster.AffinityPropagation):

# +
# Clustering using affinity propogation
# pao1_clustering = AffinityPropagation().fit(pao1_corr)
# pa14_clustering = AffinityPropagation().fit(pa14_corr)
# -

# ## Membership assignments

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
# -

# **Final method:**
# We will use <Method> because ...
#
#     Thoughts on different methods

# +
# Save membership dataframe
# pao1_membership_df.to_csv(pao1_membership_filename, sep="\t")
# pa14_membership_df.to_csv(pa14_membership_filename, sep="\t")
