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

# %load_ext autoreload
# %autoreload 2
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from core_acc_modules import paths

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

        # Distribution plot for core genes
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

# ## Examine composition of modules
#
# We expect that genes within the same operon or regulon will cluster together (i.e. be within the same module). To test this we will compare the distribution of the number of modules that contain genes within the same regulon vs the number of modules that contain random genes

# +
# Load PAO1 regulon and operon file
pao1_regulon_filename = paths.PAO1_REGULON
pao1_operon_filename = paths.PAO1_OPERON

# Load membership for specific threshold
corr_threshold = 0.9
pao1_membership_filename = f"pao1_membership_{corr_threshold}.tsv"

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", header=0, index_col=0)
# -

pao1_membership.head()

# Parse file
pao1_operon = pd.read_csv(pao1_operon_filename, index_col=0, header=0)
print(pao1_operon.shape)
pao1_operon.head()

# For each regulon/operson, select a random set of genes that are the same size at the regulon/operon
pao1_operon["Random Genes"] = pao1_operon["Length"].apply(
    lambda num_genes: pao1_membership.sample(num_genes).index.values
)

pao1_operon.head()

pao1_operon["Genes"] = pao1_operon["Genes"].str.split(";")

pao1_operon.head()

# Need to deal with genes not in membership list
pao1_operon["Genes"].apply(
    lambda list_genes: pao1_membership.loc[list_genes]["module id"].nunique()
)

pao1_membership.loc[pao1_operon.loc[12029, "Genes"]]["module id"].nunique()

pao1_membership.shape
"PA2319" in list(pao1_membership.index)
# pao1_membership['PA2319']

# +
# For each regulon/operon get the number of modules that regulon/operon genes are found in, number of modules
# that random genes are found in

# +
# Compare distributions using t-test
