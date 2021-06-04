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

# # Validate array vs RNA-seq modules
#
# The goal of this notebook is to compare the modules found using array data vs RNA-seq data. We would expect the modules to be similar

# +
# %load_ext autoreload
# %autoreload 2
import os
import scipy
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import DBSCAN, AgglomerativeClustering, AffinityPropagation
from core_acc_modules import utils, paths

np.random.seed(1)

# +
# User params
num_singular_values_log = 100

# Clustering method
# Choices: {"dbscan", "hierarchal", "affinity"}
cluster_method = "affinity"

# DBSCAN params
density_threshold = 8

# Hierarchical clustering params
hier_threshold = 8
link_dist = "average"

# Affinity params
affinity_damping = 0.6
# -

# Load expression data
pao1_array_compendium_filename = paths.ARRAY_COMPENDIUM_TO_COMPARE
pao1_rnaseq_compendium_filename = paths.RNASEQ_COMPENDIUM_TO_COMPARE

pao1_array_compendium = pd.read_csv(
    pao1_array_compendium_filename, sep="\t", header=0, index_col=0
)
pao1_rnaseq_compendium = pd.read_csv(
    pao1_rnaseq_compendium_filename, sep="\t", header=0, index_col=0
)

print(pao1_array_compendium.shape)
pao1_array_compendium.head()

print(pao1_rnaseq_compendium.shape)
pao1_rnaseq_compendium.head()

# ## Get correlation matrices

# Correlation
pao1_array_corr_original = pao1_array_compendium.corr()

# Note: Below we plotted the heatmap of the array data to confirm that it has the same issue as the RNA-seq data - there is one large cluster

# %%time
# Plot heatmap
"""o1 = sns.clustermap(pao1_array_corr_original.abs(), cmap="viridis", figsize=(20, 20))
o1.fig.suptitle("Correlation of raw PAO1 genes (array compendium)", y=1.05)"""

# Transpose compendia to be gene x sample
# Here we're interested in how genes cluster
pao1_array_compendium_T = pao1_array_compendium.T
pao1_rnaseq_compendium_T = pao1_rnaseq_compendium.T

# log transform data
pao1_array_compendium_log10 = np.log10(1 + pao1_array_compendium_T)
pao1_rnaseq_compendium_log10 = np.log10(1 + pao1_rnaseq_compendium_T)

# Apply SVD
array_U, array_s, array_Vh = np.linalg.svd(
    pao1_array_compendium_log10, full_matrices=False
)
rnaseq_U, rnaseq_s, rnaseq_Vh = np.linalg.svd(
    pao1_rnaseq_compendium_log10, full_matrices=False
)

# Convert ndarray to df to use corr()
array_U_df = pd.DataFrame(data=array_U, index=pao1_array_compendium_T.index)
rnaseq_U_df = pd.DataFrame(data=rnaseq_U, index=pao1_rnaseq_compendium_T.index)

# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
pao1_array_corr_log_spell = array_U_df.iloc[:, :num_singular_values_log].T.corr()
pao1_rnaseq_corr_log_spell = rnaseq_U_df.iloc[:, :num_singular_values_log].T.corr()

# Note: Here we plot the heatmaps to verify that the correlation of log + SPELL transformed data looks as expected (i.e. there is not a single large cluster)

# Plot heatmap
"""h1a = sns.clustermap(pao1_array_corr_log_spell.abs(), cmap="viridis", figsize=(20, 20))
h1a.fig.suptitle(
    f"log transform + SPELL corrected using {num_singular_values_log} vectors (PAO1 array)",
    y=1.05,
)"""

# Plot heatmap
"""h1b = sns.clustermap(pao1_rnaseq_corr_log_spell.abs(), cmap="viridis", figsize=(20, 20))
h1b.fig.suptitle(
    f"log transform + SPELL corrected using {num_singular_values_log} vectors (PAO1 rnaseq)",
    y=1.05,
)"""

# ## Clustering and get module membership

# Clustering using DBSCAN
if cluster_method == "dbscan":
    pao1_array_clustering = DBSCAN(eps=density_threshold).fit(pao1_array_corr_log_spell)
    pao1_rnaseq_clustering = DBSCAN(eps=density_threshold).fit(
        pao1_rnaseq_corr_log_spell
    )

# Clustering using hierarchal clustering
if cluster_method == "hierarchal":
    pao1_array_clustering = AgglomerativeClustering(
        n_clusters=None, distance_threshold=hier_threshold, linkage=link_dist
    ).fit(ppao1_array_corr_log_spell)
    pao1_rnaseq_clustering = AgglomerativeClustering(
        n_clusters=None, distance_threshold=hier_threshold, linkage=link_dist
    ).fit(pao1_rnaseq_corr_log_spell)

# Clustering using affinity propogation
if cluster_method == "affinity":
    pao1_array_clustering = AffinityPropagation(random_state=0).fit(
        pao1_array_corr_log_spell
    )
    pao1_rnaseq_clustering = AffinityPropagation(
        random_state=1, damping=affinity_damping
    ).fit(pao1_rnaseq_corr_log_spell)

# +
# Get module membership for a single threshold
# Format and save output to have columns: gene_id | group_id
pao1_array_membership_df = pd.DataFrame(
    data={"module id": pao1_array_clustering.labels_},
    index=pao1_array_corr_log_spell.index,
)

pao1_array_membership_df["module id"].value_counts()
# -

pao1_array_membership_df.head()

# +
pao1_rnaseq_membership_df = pd.DataFrame(
    data={"module id": pao1_rnaseq_clustering.labels_},
    index=pao1_rnaseq_corr_log_spell.index,
)

pao1_rnaseq_membership_df["module id"].value_counts()
# -

pao1_rnaseq_membership_df.head()

# ## Compare composition of modules
#
# For a given array module, are the genes within 1 module in the RNA-seq compendium?

# +
array_module_mapping = {}
for grp_id, grp in pao1_array_membership_df.groupby("module id"):
    grp_mapped = pao1_rnaseq_membership_df.loc[grp.index]
    array_module_mapping[grp_id] = list(grp_mapped["module id"].unique())

array_module_mapping

# +
# Are there any array modules that map to a single rnaseq module?
consistent_array_modules = []
for array_module_id, list_rnaseq_module_ids in array_module_mapping.items():
    if len(list_rnaseq_module_ids) == 1:
        print(array_module_id, list_rnaseq_module_ids)
        consistent_array_modules.append(array_module_id)

print(len(consistent_array_modules) / len(array_module_mapping))

# +
rnaseq_module_mapping = {}
for grp_id, grp in pao1_rnaseq_membership_df.groupby("module id"):
    grp_mapped = pao1_array_membership_df.loc[grp.index]
    rnaseq_module_mapping[grp_id] = list(grp_mapped["module id"].unique())

rnaseq_module_mapping

# +
# Are there any rnaseq modules that map to a single array module?
consistent_rnaseq_modules = []
for rnaseq_module_id, list_array_module_ids in rnaseq_module_mapping.items():
    if len(list_array_module_ids) == 1:
        print(rnaseq_module_id, list_array_module_ids)
        consistent_rnaseq_modules.append(rnaseq_module_id)

print(len(consistent_rnaseq_modules) / len(rnaseq_module_mapping))
# -

# **Observation:**
#
# * Only ~3% of modules are consistent using RNA-seq and array data. We would have expected more of an overlap. What is causing this?
