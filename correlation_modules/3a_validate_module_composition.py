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
# This notebook examines the composition of the modules is a few ways:
# 1. This notebook compares the modules found using array data vs RNA-seq data. We would expect the modules to be similar
# 2. This notebook also examines the enrichment of KEGG pathways within the modules

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
import matplotlib.pyplot as plt
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

pao1_array_membership_df["module id"].nunique()

# Wondering if the size of the modules is affecting the p-values
# Let's try removing the largest size modules
module_ids = pao1_array_membership_df["module id"].value_counts()
module_ids_toexclude = module_ids[module_ids > 20].index

# +
# pao1_array_membership_df["module id"] in module_ids_toexclude
# -

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
#
# To compare the consistency we use the [hypergeometric test](https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458) which will allow us to examine the over-representation of genes from the RNA-seq module within an array module and vice versa.

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
    total_array_genes = array_membership_df.shape[0]

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
                total_array_genes,
                num_rnaseq_genes,
                num_array_genes,
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
    total_rnaseq_genes = rnaseq_membership_df.shape[0]

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
                total_rnaseq_genes,
                num_array_genes,
                num_rnaseq_genes,
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

# +
# Test: mean p-values using the true module labels vs shuffle module labels is significant
true_array_to_rnaseq = array_to_rnaseq_combined[
    array_to_rnaseq_combined["label"] == "true"
]["p-value"].values
shuffled_array_to_rnaseq = array_to_rnaseq_combined[
    array_to_rnaseq_combined["label"] == "shuffle"
]["p-value"].values

(stats, pvalue) = scipy.stats.ttest_ind(true_array_to_rnaseq, shuffled_array_to_rnaseq)
print(pvalue)
# -

f = sns.boxplot(
    x=array_to_rnaseq_combined["label"],
    y=array_to_rnaseq_combined["p-value"],
    palette="Blues",
    notch=True,
)
f.set_title("Array to RNA-seq modules")

# +
# Test: mean p-values using the true module labels vs shuffle module labels is significant
true_rnaseq_to_array = rnaseq_to_array_combined[
    rnaseq_to_array_combined["label"] == "true"
]["p-value"].values
shuffled_rnaseq_to_array = rnaseq_to_array_combined[
    rnaseq_to_array_combined["label"] == "shuffle"
]["p-value"].values

(stats, pvalue) = scipy.stats.ttest_ind(true_rnaseq_to_array, shuffled_rnaseq_to_array)
print(pvalue)
# -

g = sns.boxplot(
    x=rnaseq_to_array_combined["label"],
    y=rnaseq_to_array_combined["p-value"],
    palette="Blues",
    notch=True,
)
g.set_title("RNA-seq to array modules")

# **Takeaway:**
#
# * Comparing how RNA-seq modules map to array modules, looks like most modules have very low p-values, meaning they map well between RNA-seq and array.
# * In contrast, random modules have slightly higher p-values, indicating that the modules don't map as well, though the p-values are still very low.
# * Using t-test, there is a signficant difference in the distribution of p-values using the true module labels versus the shuffled module labels.

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
rnaseq_kegg_true.head()

# %%time
rnaseq_kegg_shuffle = map_rnaseq_to_KEGG(
    pao1_rnaseq_membership_shuffle_df, pao1_pathways
)
rnaseq_kegg_shuffle.head()

rnaseq_kegg_true.describe()

rnaseq_kegg_shuffle.describe()

# +
rnaseq_kegg_true["label"] = "true"
rnaseq_kegg_shuffle["label"] = "shuffle"

rnaseq_kegg_combined = pd.concat([rnaseq_kegg_true, rnaseq_kegg_shuffle])

# +
# Test: mean p-values using the true module labels vs shuffle module labels is significant
true_rnaseq_to_kegg = rnaseq_kegg_combined[rnaseq_kegg_combined["label"] == "true"][
    "p-value"
].values
shuffled_rnaseq_to_kegg = rnaseq_kegg_combined[
    rnaseq_kegg_combined["label"] == "shuffle"
]["p-value"].values

(stats, pvalue) = scipy.stats.ttest_ind(true_rnaseq_to_kegg, shuffled_rnaseq_to_kegg)
print(pvalue)
# -

h = sns.boxplot(
    x=rnaseq_kegg_combined["label"],
    y=rnaseq_kegg_combined["p-value"],
    palette="Blues",
    notch=True,
)
h.set_title("KEGG enrichment in RNA-seq modules")

# **Takeaway:**
#
# * The p-values that correspond to the over-representation of KEGG pathways within modules is low using both the true module labels and the shuffled module labels. But the shuffled p-values are lower in this case, I'm not sure why.
