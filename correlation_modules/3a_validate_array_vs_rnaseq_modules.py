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
# %matplotlib inline
import os
import scipy
import scipy.stats as ss
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

# Load array expression data
pao1_array_compendium_filename = paths.ARRAY_COMPENDIUM_TO_COMPARE

pao1_array_compendium = pd.read_csv(
    pao1_array_compendium_filename, sep="\t", header=0, index_col=0
)

print(pao1_array_compendium.shape)
pao1_array_compendium.head()

# ## Get correlation matrix for array compendium

# Correlation
pao1_array_corr_original = pao1_array_compendium.corr()

# Note: Below we plotted the heatmap of the array data to confirm that it has the same issue as the RNA-seq data - there is one large cluster

# %%time
# Plot heatmap
o1 = sns.clustermap(pao1_array_corr_original.abs(), cmap="viridis", figsize=(20, 20))
o1.fig.suptitle("Correlation of raw PAO1 genes (array compendium)", y=1.05)

# Transpose compendia to be gene x sample
# Here we're interested in how genes cluster
pao1_array_compendium_T = pao1_array_compendium.T

# log transform data
pao1_array_compendium_log10 = np.log10(1 + pao1_array_compendium_T)

# Apply SVD
array_U, array_s, array_Vh = np.linalg.svd(
    pao1_array_compendium_log10, full_matrices=False
)

# Convert ndarray to df to use corr()
array_U_df = pd.DataFrame(data=array_U, index=pao1_array_compendium_T.index)

# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
pao1_array_corr_log_spell = array_U_df.iloc[:, :num_singular_values_log].T.corr()

# Note: Here we plot the heatmaps to verify that the correlation of log + SPELL transformed data looks as expected (i.e. there is not a single large cluster)

# %%time
# Plot heatmap
h1a = sns.clustermap(pao1_array_corr_log_spell.abs(), cmap="viridis", figsize=(20, 20))
h1a.fig.suptitle(
    f"log transform + SPELL corrected using {num_singular_values_log} vectors (PAO1 array)",
    y=1.05,
)

# ## Clustering and get module membership for array compendium

# Clustering using DBSCAN
if cluster_method == "dbscan":
    pao1_array_clustering = DBSCAN(eps=density_threshold).fit(pao1_array_corr_log_spell)

# Clustering using hierarchal clustering
if cluster_method == "hierarchal":
    pao1_array_clustering = AgglomerativeClustering(
        n_clusters=None, distance_threshold=hier_threshold, linkage=link_dist
    ).fit(ppao1_array_corr_log_spell)

# Clustering using affinity propogation
if cluster_method == "affinity":
    pao1_array_clustering = AffinityPropagation(random_state=0).fit(
        pao1_array_corr_log_spell
    )

# Get module membership for a single threshold
# Format and save output to have columns: gene_id | group_id
pao1_array_membership_df = pd.DataFrame(
    data={"module id": pao1_array_clustering.labels_},
    index=pao1_array_corr_log_spell.index,
)
pao1_array_membership_df["module id"].value_counts()

pao1_array_membership_df.head()

# ### Load RNA-seq module membership
#
# The modules for the RNA-seq compendium were generated using the same procedure as we just performed using the array compendium. The code for performing this module detection can be found in the previous 1_ and 2_ notebooks in this directory.

pao1_rnaseq_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{cluster_method}_all.tsv"
)

pao1_rnaseq_membership_df = pd.read_csv(
    pao1_rnaseq_membership_filename, sep="\t", index_col=0, header=0
)

pao1_rnaseq_membership_df.head()

# ## Compare composition of modules between array and RNA-seq
#
# For a given RNA-seq module, find the best fit module in the array compendium using the hypergeometric test?
#
# For rna-seq data, `modules id` is the cluster id. We can then look up which array cluster these genes from the rna-seq cluster map to. Here we are looking to see if modules are consistent (i.e. Do genes that cluster together in the RNA-seq compendium also cluster in the array compendium)?

# +
# As a baseline make a membership df mapping genes to a shuffled set of module ids
pao1_rnaseq_membership_shuffle_df = pao1_rnaseq_membership_df.copy()
pao1_rnaseq_membership_shuffle_df["module id"] = np.random.permutation(
    pao1_rnaseq_membership_shuffle_df["module id"].values
)

pao1_array_membership_shuffle_df = pao1_array_membership_df.copy()
pao1_array_membership_shuffle_df["module id"] = np.random.permutation(
    pao1_array_membership_shuffle_df["module id"].values
)
# -

pao1_rnaseq_membership_shuffle_df.head()

pao1_array_membership_shuffle_df.head()


# Array to RNA-seq
# Given an array module, look for the rnaseq module with the most overlap/significant p-value
def map_array_to_rnaseq_modules(array_membership_df, rnaseq_membership_df):
    total_rnaseq_genes = rnaseq_membership_df.shape[0]

    rows = []
    # For each array module
    for array_group_name, array_df_group in array_membership_df.groupby("module id"):
        num_array_genes = array_df_group.shape[0]

        # Find the RNA-seq module with the best overlap
        rnaseq_mapping = {}
        for rnaseq_group_name, rnaseq_df_group in rnaseq_membership_df.groupby(
            "module id"
        ):
            num_rnaseq_genes = rnaseq_df_group.shape[0]

            shared_genes = set(array_df_group.index).intersection(rnaseq_df_group.index)
            num_shared_genes = len(shared_genes)

            pval = ss.hypergeom.sf(
                num_shared_genes,
                total_rnaseq_genes,
                num_array_genes,
                num_rnaseq_genes,
            )

            # rnaseq_mapping[rnaseq_module] = p-value
            rnaseq_mapping[rnaseq_group_name] = pval

        # Find the module with the best mapping (i.e. lowest p-value)
        best_rnaseq_module = min(rnaseq_mapping, key=rnaseq_mapping.get)

        rows.append(
            {
                "array module id": array_group_name,
                "best rnaseq module id": best_rnaseq_module,
                "p-value": min(rnaseq_mapping.values()),
            }
        )

    array_to_rnaseq_df = pd.DataFrame(rows)

    return array_to_rnaseq_df


# %%time
array_to_rnaseq_true = map_array_to_rnaseq_modules(
    pao1_array_membership_df, pao1_rnaseq_membership_df
)

array_to_rnaseq_true.head()

# %%time
array_to_rnaseq_shuffle = map_array_to_rnaseq_modules(
    pao1_array_membership_df, pao1_rnaseq_membership_shuffle_df
)

array_to_rnaseq_shuffle.head()


# RNA-seq to array
# Given a RNA-seq module, look for the array module with the most overlap/significant p-value
def map_rnaseq_to_array_modules(rnaseq_membership_df, array_membership_df):
    total_array_genes = array_membership_df.shape[0]

    rows = []
    # For each rna-seq module
    for rnaseq_group_name, rnaseq_df_group in rnaseq_membership_df.groupby("module id"):
        num_rnaseq_genes = rnaseq_df_group.shape[0]

        # Find the array module with the best overlap
        array_mapping = {}
        for array_group_name, array_df_group in array_membership_df.groupby(
            "module id"
        ):
            num_array_genes = array_df_group.shape[0]

            shared_genes = set(array_df_group.index).intersection(rnaseq_df_group.index)
            num_shared_genes = len(shared_genes)

            pval = ss.hypergeom.sf(
                num_shared_genes,
                total_array_genes,
                num_rnaseq_genes,
                num_array_genes,
            )

            # array_mapping[rnaseq_module] = p-value
            array_mapping[array_group_name] = pval

        # Find the module with the best mapping (i.e. lowest p-value)
        best_array_module = min(array_mapping, key=array_mapping.get)

        rows.append(
            {
                "RNA-seq module id": rnaseq_group_name,
                "best array module id": best_array_module,
                "p-value": min(array_mapping.values()),
            }
        )

    rnaseq_to_array_df = pd.DataFrame(rows)

    return rnaseq_to_array_df


# %%time
rnaseq_to_array_true = map_rnaseq_to_array_modules(
    pao1_rnaseq_membership_df, pao1_array_membership_df
)

rnaseq_to_array_true.head()

# %%time
rnaseq_to_array_shuffle = map_rnaseq_to_array_modules(
    pao1_rnaseq_membership_df, pao1_array_membership_shuffle_df
)

rnaseq_to_array_shuffle.head()

# ## Plot

# +
# Format df for plotting
array_to_rnaseq_true["label"] = "true"
array_to_rnaseq_shuffle["label"] = "shuffle"

array_to_rnaseq_combined = pd.concat([array_to_rnaseq_true, array_to_rnaseq_shuffle])

# +
rnaseq_to_array_true["label"] = "true"
rnaseq_to_array_shuffle["label"] = "shuffle"

rnaseq_to_array_combined = pd.concat([rnaseq_to_array_true, rnaseq_to_array_shuffle])
# -

array_to_rnaseq_combined.head()

f = sns.boxplot(
    x=array_to_rnaseq_combined["label"],
    y=array_to_rnaseq_combined["p-value"],
    palette="Blues",
)
f.set_title("Array to RNA-seq modules")

g = sns.boxplot(
    x=rnaseq_to_array_combined["label"],
    y=rnaseq_to_array_combined["p-value"],
    palette="Blues",
)
g.set_title("RNA-seq to array modules")

# **Observation:**
#
# * Comparing how RNA-seq modules map to array modules, looks like most modules have very low p-values, meaning they map well.
# * In contrast, random modules have slightly higher p-values, indicating that the modules don't map as well.

# ## Enrichment of modules in KEGG pathways

# +
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", index_col=0, header=None)
# -

pao1_pathways[2] = pao1_pathways[2].str.split(";").apply(set)
pao1_pathways.index = pao1_pathways.index.str.split(" - ").str[0]
pao1_pathways.head()


# RNA-seq to KEGG pathways
# Given a RNA-seq module, look for the array module with the most overlap/significant p-value
def map_rnaseq_to_KEGG(rnaseq_membership_df, kegg_df):
    total_rnaseq_genes = rnaseq_membership_df.shape[0]
    rows = []
    # For each rna-seq module
    for rnaseq_group_name, rnaseq_df_group in rnaseq_membership_df.groupby("module id"):
        num_rnaseq_genes = rnaseq_df_group.shape[0]

        # Find the KEGG pathway with the best overlap
        kegg_mapping = {}
        for kegg_name in kegg_df.index:
            num_kegg_genes = kegg_df.loc[kegg_name, 1]
            kegg_genes = list(kegg_df.loc[kegg_name, 2])

            shared_genes = set(rnaseq_df_group.index).intersection(kegg_genes)
            num_shared_genes = len(shared_genes)

            pval = ss.hypergeom.sf(
                num_shared_genes,
                total_rnaseq_genes,
                num_kegg_genes,
                num_rnaseq_genes,
            )

            # array_mapping[rnaseq_module] = p-value
            kegg_mapping[kegg_name] = pval

        # Find the module with the best mapping (i.e. lowest p-value)
        best_kegg_pathway = min(kegg_mapping, key=kegg_mapping.get)

        rows.append(
            {
                "RNA-seq module id": rnaseq_group_name,
                "best array module id": best_kegg_pathway,
                "p-value": min(kegg_mapping.values()),
            }
        )

    rnaseq_to_kegg_df = pd.DataFrame(rows)

    return rnaseq_to_kegg_df


# %%time
rnaseq_kegg_true = map_rnaseq_to_KEGG(pao1_rnaseq_membership_df, pao1_pathways)

# %%time
rnaseq_kegg_shuffle = map_rnaseq_to_KEGG(
    pao1_rnaseq_membership_shuffle_df, pao1_pathways
)

# +
rnaseq_kegg_true["label"] = "true"
rnaseq_kegg_shuffle["label"] = "shuffle"

rnaseq_kegg_combined = pd.concat([rnaseq_kegg_true, rnaseq_kegg_shuffle])
# -

h = sns.boxplot(
    x=rnaseq_kegg_combined["label"],
    y=rnaseq_kegg_combined["p-value"],
    palette="Blues",
)
h.set_title("KEGG enrichment in RNA-seq modules")

# https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
#
# Add description of hypergeometric test
