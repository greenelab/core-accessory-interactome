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

# # Explore exoS/exoU modules
#
# This notebook specifically explores the exoS (PAO1) and exoU (PA14) accessory-accessory modules to determine if there is an interesting biological story here.
#
# _P. aeruginosa_ uses a type III secretion system (T3SS) to promote development of severe disease, particularly in patients with impaired immune defenses. _P. aeruginosa_ uses a type III secretion system to inject toxic effector proteins into the cytoplasm of eukaryotic cells. ExoU, ExoS, and ExoT, three effector proteins secreted by this system. ExoU and ExoS are usually secreted by different strains.
#
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC529154/

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import utils, paths

np.random.seed(1)

# +
# Clustering method
method_name = "affinity"

# Gene subset
gene_subset = "acc"

# Select modules containing exoS (module 16) and exoU (module 17)
exoS_module_id = 16
exoU_module_id = 17
# -

# ### Load correlation matrix

# Load correlation matrix
pao1_corr_filename = paths.PAO1_CORR_LOG_SPELL_ACC
pa14_corr_filename = paths.PA14_CORR_LOG_SPELL_ACC

pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pa14_corr = pd.read_csv(pa14_corr_filename, sep="\t", index_col=0, header=0)

# ### Load module membership

pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method_name}_{gene_subset}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{method_name}_{gene_subset}.tsv"
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", header=0, index_col=0)
pa14_membership = pd.read_csv(pa14_membership_filename, sep="\t", header=0, index_col=0)

# ### Select genes associated with modules of interest

exoS_module_df = pao1_membership[pao1_membership["module id"] == exoS_module_id]
exoU_module_df = pa14_membership[pa14_membership["module id"] == exoU_module_id]

exoS_module_df.head()

exoU_module_df.head()

# ### Get most co-expressed accessory genes
#
# Make dataframe with all accessory genes sorted by their correlation score relative to exoS or exoU

# Get gene id for exoS and exoU
exoS_id = "PA3841"
exoU_id = "PA14_51530"

# Get correlation scores
exoS_corr_all = pao1_corr.loc[exoS_id].to_frame("corr_score")
exoU_corr_all = pa14_corr.loc[exoU_id].to_frame("corr_score")

exoS_corr_all.head()

exoU_corr_all.head()

# +
# Import gene metadata
pao1_gene_annot_filename = paths.GENE_PAO1_ANNOT
pa14_gene_annot_filename = paths.GENE_PA14_ANNOT

pao1_gene_annot = pd.read_csv(pao1_gene_annot_filename, index_col=0, header=0)
pa14_gene_annot = pd.read_csv(pa14_gene_annot_filename, index_col=0, header=0)
# -

# Get df with gene ids as indices and gene names as a column
# Having the data in a df instead of a series will just allow me to do my merges that are in the notebook
pao1_gene_annot = pao1_gene_annot["Name"].to_frame("gene name")
pa14_gene_annot = pa14_gene_annot["Name"].to_frame("gene name")

# Make dataframe with gene id, correlation score, gene name
# Add gene names
exoS_corr_df = exoS_corr_all.merge(
    pao1_gene_annot, left_index=True, right_index=True, how="left"
)
exoU_corr_df = exoU_corr_all.merge(
    pa14_gene_annot, left_index=True, right_index=True, how="left"
)

print(exoS_corr_df.shape)
exoS_corr_df.head()

print(exoU_corr_df.shape)
exoU_corr_df.head()

# Save
exoS_corr_df.to_csv("exoS_corr.tsv", sep="\t")
exoU_corr_df.to_csv("exoU_corr.tsv", sep="\t")

# ### Heatmap

exoS_module_gene_ids = list(exoS_module_df.index)
exoU_module_gene_ids = list(exoU_module_df.index)

exoS_corr = pao1_corr.loc[exoS_module_gene_ids, exoS_module_gene_ids]
exoU_corr = pa14_corr.loc[exoU_module_gene_ids, exoU_module_gene_ids]

# %%time
f = sns.clustermap(exoS_corr.abs(), cmap="viridis", figsize=(20, 20))
f.ax_heatmap.set_xticklabels(f.ax_heatmap.get_xmajorticklabels(), fontsize=20)
f.ax_heatmap.set_yticklabels(f.ax_heatmap.get_ymajorticklabels(), fontsize=20)
f.fig.suptitle("Correlation of exoS module", y=1.05, fontsize=24)

# %%time
g = sns.clustermap(exoU_corr.abs(), cmap="viridis", figsize=(20, 20))
g.ax_heatmap.set_xticklabels(f.ax_heatmap.get_xmajorticklabels(), fontsize=20)
g.ax_heatmap.set_yticklabels(f.ax_heatmap.get_ymajorticklabels(), fontsize=20)
g.fig.suptitle("Correlation of exoU module", y=1.05, fontsize=24)

# **Takeaway**
# We've aggregated the information from this notebook into a [google sheet](https://docs.google.com/spreadsheets/d/1AuD1Q4lHhWNp5xzgW-hi8mHkHFyd91rmOksXXuAwk4Q/edit#gid=533448426) to easily share with collaborators. This sheet contains Uniprot annotations for each gene within the exoS and exoU modules. The sheet also contains a sorted matrix of genes and how correlated they are with exoS and exoU.
#
# * Genes within exoS module appear to be more highly correlated with exoS (see sheet) and each other (see heatmap) unlike exoU module
#     * What might this mean about exoS, exoU?
#     * Despite being part of the same T3S system, their relationship to other accessory genes is different
# * Based on gene annotations, is there a different mechanism by which exoS contributes to virulence compared to exoU?
#     * This is difficult to answer with so many unannotated genes
#
# Some more reading will need to be done to determine the biological motivation here:
# * What is known about the mechanism by which these genes contribute to virulence?
# * What can we learn from module composition?
# * What can we learn from most co-expressed genes?
