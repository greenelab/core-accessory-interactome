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

# # Correlation analysis
#
# This notebook creates the correlation matrix using the MR counts and a SPELL processed version.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import plotnine as pn
import seaborn as sns
from sklearn import preprocessing
import matplotlib.pyplot as plt
import umap
import random
import numpy as np
import scipy
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from scripts import paths, utils

# ## Set user parameters
#
# Here we set the number of PCs or singular vectors to use. We are starting with 300 since this is what [eADAGE](https://pubmed.ncbi.nlm.nih.gov/28711280/) used.

# +
# Params

# Which subset of genes to consider: core, acc, all
subset_genes = "all"

# The number of accessory genes is 200 - 500
# The number of core genes is ~ 5000
# So the number of singular vectors is relative to the number of genes
# These numbers for selected based on a manual inspection of the heatmaps
# Making sure the dominant signal is removed and the approximate size of the clusters
# seems reasonable
if subset_genes == "acc":
    num_SVs = 50
else:
    num_SVs = 100
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

# ## Get core/accessory genes
#
# We will subset the correlation matrix to only consider core genes.
#
# _Rationale:_ Previously we used all genes (both core and accessory) to create a co-expression network, but due to the large imbalance in the number of core genes compared to accessory genes, no module was found to be "mostly core." Instead we will perform separate analyses of core and accessory genes to examine co-expression patterns.

# +
# Read in expression data
pao1_expression_filename = paths.PAO1_COMPENDIUM
pa14_expression_filename = paths.PA14_COMPENDIUM

pao1_expression = pd.read_csv(pao1_expression_filename, sep="\t", index_col=0, header=0)
pa14_expression = pd.read_csv(pa14_expression_filename, sep="\t", index_col=0, header=0)

# +
pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT

core_acc_dict = utils.get_my_core_acc_genes(
    pao1_annot_filename, pa14_annot_filename, pao1_expression, pa14_expression
)
# -

pao1_core = core_acc_dict["core_pao1"]
pa14_core = core_acc_dict["core_pa14"]
pao1_acc = core_acc_dict["acc_pao1"]
pa14_acc = core_acc_dict["acc_pa14"]

# ## Select gene subset

# Select subset of genes
if subset_genes == "core":
    pao1_compendium = pao1_compendium[pao1_core]
    pa14_compendium = pa14_compendium[pa14_core]
elif subset_genes == "acc":
    pao1_compendium = pao1_compendium[pao1_acc]
    pa14_compendium = pa14_compendium[pa14_acc]

print(pao1_compendium.shape)
print(pa14_compendium.shape)

# ## Correlation of raw gene expression data
#
# Here is the correlation of the raw data without any manipulations. This will serve as a reference to compare the correlations below where applied corrections to the correlations to account for the dominant signal described above.

# Correlation
pao1_corr_original = pao1_compendium.corr()
pa14_corr_original = pa14_compendium.corr()

# Check for duplicates indices
assert pao1_corr_original.index.duplicated().sum() == 0
assert pa14_corr_original.index.duplicated().sum() == 0

# Check for duplicate rows
assert pao1_corr_original[pao1_corr_original.duplicated(keep=False)].shape[0] == 0
assert pa14_corr_original[pa14_corr_original.duplicated(keep=False)].shape[0] == 0

print(pao1_corr_original.shape)
pao1_corr_original.head()

print(pa14_corr_original.shape)
pa14_corr_original.head()

# +
# %%time
# Plot heatmap
o1 = sns.clustermap(pao1_corr_original, cmap="BrBG", center=0, figsize=(20, 20))
o1.fig.suptitle("Correlation of raw PAO1 genes", y=1.05, fontsize=24)

# Save
pao1_pearson_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_raw_clustermap_test.png"
)
o1.savefig(pao1_pearson_filename, dpi=300)

# +
# Plot heatmap
o2 = sns.clustermap(pa14_corr_original, cmap="BrBG", center=0, figsize=(20, 20))
o2.fig.suptitle("Correlation of raw PA14 genes", y=1.05, fontsize=24)

# Save
pa14_pearson_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14__{subset_genes}_raw_clustermap_test.png"
)
o2.savefig(pa14_pearson_filename, dpi=300)
# -

# ## Log transform + SPELL Correlation
#
# _Review of SVD_
#
# Singular Value Decomposition is a way to factorize your matrix, $X^{mxn}$ into singular vectors and singular values: $X = U \Sigma V^*$
#
# In our case $X$ is **gene x sample** and then the columns of $U$ (gene x eigensample) are the left singular vectors (gene coefficient vectors); $\Sigma$ (eigengene x eigensample) has singular values and is diagonal (mode amplitudes); and $V^T$ (eigengene x sample) has rows that are the right singular vectors (expression level vectors).
#
# Here we are using SVD to reduce the noise in our original data by performing dimensionality reduction. This dimensionality is done by neglecting the small singular values in the diagonal matrix $\Sigma$. Normally people would get the projection of the original data onto the singular vectors by $U \Sigma$ and apply the correlation on the projected data. Here, we're following the description in [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926) where they performed correlation on $U$ only.
#
# From [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926), we apply their "signal balancing technique that enhances biological information". This is the first part of their [SPELL](https://spell.yeastgenome.org/) algorithm that is described in section 2.3.1. SPELL calculates the correlation on the gene coefficient matrix, $U$ (i.e. how much genes contribute to a latent variable) that is generated after applying SVD. This matrix represents how genes contribute to independent latent variables that capture the signal in the data where the variance of the variables is 1. The idea is that correlations between gene contributions are more balanced so that less prominent patterns are amplified and more dominant patterns are dampended due to this compression. Figure 3 shows how well SPELL recapitulates biology (i.e. the relationship between genes within a GO term) compared to Pearson correlation.

# Transpose compendia to be gene x sample
# Here we're interested in how genes cluster
pao1_compendium_T = pao1_compendium.T
pa14_compendium_T = pa14_compendium.T

# log transform data
pao1_compendium_log10 = np.log10(1 + pao1_compendium_T)
pa14_compendium_log10 = np.log10(1 + pa14_compendium_T)

# Apply SVD
pao1_U, pao1_s, pao1_Vh = np.linalg.svd(pao1_compendium_log10, full_matrices=False)
pa14_U, pa14_s, pa14_Vh = np.linalg.svd(pa14_compendium_log10, full_matrices=False)

# #### Quick check
#
# Plot the variance explained to make sure that our choice of number of singular vectors is reasonable and aligns with our manual inspection

plt.plot(pao1_s ** 2 / sum(pao1_s ** 2) * 100)
plt.ylabel("Percent variability explained")

plt.plot(pa14_s ** 2 / sum(pa14_s ** 2) * 100)
plt.ylabel("Percent variability explained")

print(pao1_compendium_T.shape)
print(pao1_U.shape, pao1_s.shape, pao1_Vh.shape)

print(pa14_compendium_T.shape)
print(pa14_U.shape, pa14_s.shape, pa14_Vh.shape)

# Convert ndarray to df to use corr()
pao1_U_df = pd.DataFrame(data=pao1_U, index=pao1_compendium_T.index)
pa14_U_df = pd.DataFrame(data=pa14_U, index=pa14_compendium_T.index)

# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
pao1_corr_log_spell = pao1_U_df.iloc[:, :num_SVs].T.corr()
pa14_corr_log_spell = pa14_U_df.iloc[:, :num_SVs].T.corr()

# Check for duplicates indices
assert pao1_corr_log_spell.index.duplicated().sum() == 0
assert pa14_corr_log_spell.index.duplicated().sum() == 0

# Check for duplicate rows
assert pao1_corr_log_spell[pao1_corr_log_spell.duplicated(keep=False)].shape[0] == 0
assert pa14_corr_log_spell[pa14_corr_log_spell.duplicated(keep=False)].shape[0] == 0

# +
# Plot heatmap
h1a = sns.clustermap(pao1_corr_log_spell, cmap="BrBG", center=0, figsize=(20, 20))
h1a.fig.suptitle(
    f"log transform + SPELL corrected using {num_SVs} vectors (PAO1)",
    y=1.05,
    fontsize=24,
)

# Save
pao1_log_spell_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_log_spell_clustermap_test.png"
)
h1a.savefig(pao1_log_spell_filename, dpi=300)

# +
h2a = sns.clustermap(pa14_corr_log_spell, cmap="BrBG", center=0, figsize=(20, 20))
h2a.fig.suptitle(
    f"log transformed + SPELL corrected using {num_SVs} vectors (PA14)",
    y=1.05,
    fontsize=24,
)

# Save
pa14_log_spell_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_{subset_genes}_log_spell_clustermap_test.png"
)
h2a.savefig(pa14_log_spell_filename, dpi=300)
# -

# Save raw correlation matrices
pao1_original_mat_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_raw_mat_test.tsv"
)
pa14_original_mat_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_{subset_genes}_raw_mat_test.tsv"
)
pao1_corr_original.to_csv(pao1_original_mat_filename, sep="\t")
pa14_corr_original.to_csv(pa14_original_mat_filename, sep="\t")

# Save log transform + SPELL correlation matrices
pao1_log_spell_mat_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_log_spell_mat_test.tsv"
)
pa14_log_spell_mat_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_{subset_genes}_log_spell_mat_test.tsv"
)
pao1_corr_log_spell.to_csv(pao1_log_spell_mat_filename, sep="\t")
pa14_corr_log_spell.to_csv(pa14_log_spell_mat_filename, sep="\t")
