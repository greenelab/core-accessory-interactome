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

# ## Quick validation of network modules

# +
# %load_ext autoreload
# %autoreload 2
import os
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from core_acc_modules import paths

np.random.seed(1)

# +
# User params
# Remove modules of this size or greater for analysis looking at coverage of regulon/operons
module_size_threshold = 1000

# Correlation threshold to determine which network modules to examine
corr_threshold_toexamine = 0.9
# -

# ## Examine size of modules

# +
corr_threshold_list = [0.9, 0.8, 0.7, 0.6, 0.5]

for corr_threshold in corr_threshold_list:
    print(f"Modules using correlation threshold: {corr_threshold}")
    pao1_membership_filename = f"pao1_membership_{corr_threshold}.tsv"
    pa14_membership_filename = f"pa14_membership_{corr_threshold}.tsv"

    pao1_membership = pd.read_csv(
        pao1_membership_filename, sep="\t", header=0, index_col=0
    )
    pa14_membership = pd.read_csv(
        pa14_membership_filename, sep="\t", header=0, index_col=0
    )

    print(pao1_membership["module id"].value_counts())
    print(pa14_membership["module id"].value_counts())


# -

# plotting function
def plot_dist_modules(threshold_list):

    # Set up the matplotlib figure
    fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(15, 15))
    axes = axes.ravel()

    for i in range(len(threshold_list)):
        pao1_membership_filename = f"pao1_membership_{threshold_list[i]}.tsv"
        pa14_membership_filename = f"pa14_membership_{threshold_list[i]}.tsv"

        pao1_membership = pd.read_csv(
            pao1_membership_filename, sep="\t", header=0, index_col=0
        )
        pa14_membership = pd.read_csv(
            pa14_membership_filename, sep="\t", header=0, index_col=0
        )

        # Get bins using all data
        hist, bins_corr = np.histogram(
            np.concatenate([pao1_membership["module id"], pa14_membership["module id"]])
        )

        fig = sns.distplot(
            pao1_membership["module id"],
            label="PAO1 modules",
            color="red",
            bins=bins_corr,
            kde=False,
            ax=axes[i],
        )

        fig = sns.distplot(
            pa14_membership["module id"],
            label="PA14 modules",
            color="blue",
            bins=bins_corr,
            kde=False,
            ax=axes[i],
        )

        fig.set_title(
            f"Histogram of size of modules using threshold {threshold_list[i]}",
            fontsize=12,
        )
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles, labels, loc="center right")


# Plot distribution of modules per threshold
plot_dist_modules(corr_threshold_list)

# **Takeaway:**
# * Looks like as we decrease our correlation threshold (i.e. connections don't need to be very strong between nodes), more modules are able to form which is what I would expect.
# * Genes in PA14 are less concentrated in one large module - there are many more small modules found.

# ## Examine composition of modules
#
# We expect that genes within the same operon or regulon will cluster together (i.e. be within the same module). To test this we will compare the distribution of the number of modules that contain genes within the same regulon vs the number of modules that contain random genes
#
# _Some definitions:_
#
# [Operons](https://en.wikipedia.org/wiki/Operon#:~:text=An%20operon%20is%20made%20up,transcription%20of%20the%20structural%20genes.) are a group of genes that share a promoter (DNA sequence that is recognized by RNA polymerase and enables transcription) and an operator (DNA sequence that repressor binds to and blocks RNA polymerase). Therefore these group of genes are transcribed or turned off together (so we would expect a very high correlation amongst these genes)
#
# [Regulons](https://en.wikipedia.org/wiki/Regulon) are a group of genes that are regulated by the same regulatory protein. A regulon can be composed of multiple operons.

# +
# Load PAO1 regulon and operon file
pao1_regulon_filename = paths.PAO1_REGULON
pao1_operon_filename = paths.PAO1_OPERON

# Load membership for specific threshold
pao1_membership_filename = f"pao1_membership_{corr_threshold_toexamine}.tsv"

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", header=0, index_col=0)
# -

pao1_membership.head()

# According to Jake relationships tend to be more meaningful if the module is smaller (e.g. if an operon with 5 genes is contained in a module consisting of 10 total genes, this seems more biologically/functionally meaningful than an operon with 5 genes contained in a module consisting of 500 genes).
#
# To correct for the single or couple very large modules, we will remove them from the analysis

# +
# Get module ids that exceed size limit
module_todrop = (
    pao1_membership["module id"]
    .value_counts()[
        (pao1_membership["module id"].value_counts() > module_size_threshold)
    ]
    .index
)

print(module_todrop)

# +
# Get genes to drop
genes_todrop = pao1_membership[pao1_membership["module id"].isin(module_todrop)].index

# Drop genes
pao1_membership = pao1_membership.drop(genes_todrop)
# -

# ### Format operon/regulon files
#
# * Remove genes from operons/regulons that don't have membership information
# * Make random list of genes with matched size
# * There are many single gene operons, we will remove these for this analysis

# +
# Read file
pao1_operon = pd.read_csv(pao1_operon_filename, index_col=0, header=0)
pao1_regulon = pd.read_csv(pao1_regulon_filename, index_col=0, header=0)

print(pao1_operon.shape)
pao1_operon.head()
# -

print(pao1_regulon.shape)
pao1_regulon.head()

# Convert "Genes" column from str to list
pao1_operon["Genes"] = pao1_operon["Genes"].str.split(";")
pao1_regulon["Genes"] = pao1_regulon["Genes"].str.split(";")

# Check if genes within operon/regulon have membership information
# Only keep genes that are found in "pao1_membership"
pao1_operon["Genes_processed"] = pao1_operon["Genes"].apply(
    lambda list_genes: [
        gene_id for gene_id in list_genes if gene_id in pao1_membership.index
    ]
)
pao1_regulon["Genes_processed"] = pao1_regulon["Genes"].apply(
    lambda list_genes: [
        gene_id for gene_id in list_genes if gene_id in pao1_membership.index
    ]
)

# Update length based on filtered gene list ("Genes_processed" column)
pao1_operon["Length_processed"] = pao1_operon["Genes_processed"].str.len()
pao1_regulon["Length_processed"] = pao1_regulon["Genes_processed"].str.len()

# +
# Quick look at distribution of size of regulons and operons
# Update length based on filtered gene list ("Genes_processed" column)

fig, axes = plt.subplots(ncols=2, nrows=1)

fig = sns.distplot(
    pao1_operon["Length_processed"],
    label="PAO1 operon size",
    color="red",
    kde=False,
    ax=axes[0],
)

fig = sns.distplot(
    pao1_regulon["Length_processed"],
    label="PAO1 regulon size",
    color="blue",
    kde=False,
    ax=axes[1],
)

fig.set_title(
    "Histogram of size of operons/regulons after filtering by membership",
    fontsize=12,
)

# +
# If number genes in operon are 1 then remove
# Drop operons and regulons that have 0 genes due to no module filtering
pao1_operon = pao1_operon.drop(pao1_operon.query("Length_processed<=1").index)
pao1_regulon = pao1_regulon.drop(pao1_regulon.query("Length_processed<=1").index)

print(pao1_operon.shape)
print(pao1_regulon.shape)
# -

# For each regulon/operon, select a random set of genes that are the same size at the regulon/operon
pao1_operon["Random_Genes"] = pao1_operon["Length_processed"].apply(
    lambda num_genes: pao1_membership.sample(num_genes).index.values
)
pao1_regulon["Random_Genes"] = pao1_regulon["Length_processed"].apply(
    lambda num_genes: pao1_membership.sample(num_genes).index.values
)

pao1_operon.head()

pao1_regulon.head()

# ### Calculate the distribution

# For each regulon/operon get the number of modules that regulon/operon genes are found in, number of modules
# that random genes are found in
pao1_operon["Num_operon_modules"] = pao1_operon["Genes_processed"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].nunique()
)
pao1_operon["Num_random_modules"] = pao1_operon["Random_Genes"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].nunique()
)

pao1_regulon["Num_regulon_modules"] = pao1_regulon["Genes_processed"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].nunique()
)
pao1_regulon["Num_random_modules"] = pao1_regulon["Random_Genes"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].nunique()
)

pao1_operon.head()

pao1_regulon.head()

# +
# Format df for plotting using displot
pao1_operon_toplot = pd.melt(
    pao1_operon, value_vars=["Num_operon_modules", "Num_random_modules"]
)
pao1_regulon_toplot = pd.melt(
    pao1_regulon, value_vars=["Num_regulon_modules", "Num_random_modules"]
)

pao1_operon_toplot.tail()


# -

def cumulative_distribution(
    data,
    scaled=False,
    survival=False,
    label="Cumulative",
    fill=False,
    flip=False,
    preserve_ends=0,
    **kwargs,
):
    """
    plots cumulative (or survival) step distribution
    adapted from https://github.com/MarvinT/morphs/blob/master/morphs/plot/utils.py
    """
    data = np.sort(data)
    if survival:
        data = data[::-1]
    y = np.arange(data.size + 1, dtype=float)
    if scaled:
        y /= y[-1]
    x = np.concatenate([data, data[[-1]]])
    plt.step(x, y, label=label, **kwargs)
    if fill:
        plt.fill_between(x, y, alpha=0.5, step="pre", **kwargs)


# +
cumulative_distribution(
    pao1_operon["Num_operon_modules"],
    label="Number of modules containing operon genes",
    color="red",
)
cumulative_distribution(
    pao1_operon["Num_random_modules"],
    label="Number of modules containing random genes",
    color="blue",
)
_ = plt.legend()
plt.title("Cumulative distribution of module counts (operon vs random genes)")
plt.ylabel("Number of operons/random groups")
plt.xlabel("The number of modules that operons/random genes are contained in")

scipy.stats.ks_2samp(
    pao1_operon["Num_operon_modules"], pao1_operon["Num_random_modules"]
)
# -

fig = sns.displot(
    pao1_operon_toplot,
    x="value",
    hue="variable",
)
plt.title("PMF distribution of module counts (operon vs random genes)")

# +
cumulative_distribution(
    pao1_regulon["Num_regulon_modules"],
    label="Number of modules containing regulon genes",
    color="red",
)
cumulative_distribution(
    pao1_regulon["Num_random_modules"],
    label="Number of modules containing random genes",
    color="blue",
)
_ = plt.legend()
plt.title("Cumulative distribution of module counts (regulon vs random genes)")
plt.ylabel("Number of regulons/random groups")
plt.xlabel("The number of modules that regulons/random genes are contained in")

scipy.stats.ks_2samp(
    pao1_regulon["Num_regulon_modules"], pao1_regulon["Num_random_modules"]
)
# -

fig = sns.displot(
    pao1_regulon_toplot,
    x="value",
    hue="variable",
)
plt.title("PMF distribution of module counts (regulon vs random genes)")

# ### Check what size the modules are that regulon/operon/random genes are found in

pao1_operon["operon_module_ids"] = pao1_operon["Genes_processed"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].unique()
)
pao1_operon["random_module_ids"] = pao1_operon["Random_Genes"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].unique()
)
pao1_operon.head()

pao1_membership["module id"].value_counts()

pao1_regulon["regulon_module_ids"] = pao1_regulon["Genes_processed"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].unique()
)
pao1_regulon["random_module_ids"] = pao1_regulon["Random_Genes"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].unique()
)
pao1_regulon.head()

# _About cumulative distribution plots:_
# * The axis cumulative distribution plots are:
#     * y-axis = The number of operon/regulon (red) or random (blue) groups.
#     * x-axis = The number of modules that operon/regulon/random genes are contained in
#
# * Looking at the operon plot, the value at "1" on the x-axis says that there are 75 operons found in exactly 1 module, and there are 10 random groups found in exactly 1 module
#     * Then the increase at "2" on the x-axis is the number of operons or random genes that are spread across 2 or 1 different modules (this is the cumulative part). In other words, a random set (size matched with the operon) where the genes in that set are found in 2 modules.
#     * There are 150 operons that are found in either 1 or 2 modules. There are 120 random groups that are found in either 1 or 2 modules.
#     * Then if you compare the blue and the red curves, the vertical distance between the two curves tells you how much of a shift there is between the distributions.
# * These distribution plots are summing counts as you move from left to right, so a shift in the curves corresponds to a shift in the distribution (i.e. a curve shifted to the right means that the distribution is shifted to the right)
#
# **Takeaway:**
# * We can perform [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) to compare the distribution of module counts for genes in regulons/operons versus random genes. The KS test will quantify the difference in the cumulative distribution curves.
#
# * Based on the KS test, there is a significant difference (across thresholds) between the operon and random distribution, as we would expect.
#     * There is only a significant difference between the regulon and random distributions at lower thresholds (0.7 and below). The lack of significance at higher thresholds is likely due to the small sample size (i.e. 6 or 10 regulons have genes that are contained in modules with fewer than 1000 genes)
