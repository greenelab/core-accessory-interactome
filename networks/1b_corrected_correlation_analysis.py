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
# This notebook examines the correlation structure in the gene expression data generated in [1_create_compendia.ipynb](../processing/1_create_compendia.ipynb).
#
# When we performed clustering on the correlation matrices (using Pearson correlation) we found that pairs of genes had either very high correlation scores (>0.5) or very low correlation scores (<0.1). As a result gene pairs that were highly correlated clustered into a single large module. This clustering pattern is not ideal for a couple of reasons:
# 1. Given that our goal is to examine the composition of gene groups, having all genes within one module does not allow us to do this
# 2. These highly correlated group of genes are likely masking other relevant specific signals/relationships in the data
#
# Here we perform two approaches to extracting correlations between genes that correct for this:
# 1. From [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926), we apply their "signal balancing technique that enhances biological information". This is the first part of their [SPELL](https://spell.yeastgenome.org/) algorithm that is described in section 2.3.1. They apply SVD to the gene expression matrix to reduce the noise that can lead to spurious results and then apply correlation between genes using the singular vector.  Correlations between genes in U equally weight each dimension of the orthonormal basis and balance their contributions such that the least prominent patterns are amplified and more dominant patterns are dampened. This process helps reveal biological signals, as some of the dominant patterns in many microarray datasets are not biologically meaningful. These dominant signals might be due to technical artifacts that lead to artificially high correlation levels (i.e. If a set of samples have all genes highly expressed or lowly expressed). Section 2.4 evaluates SPELL on its ability to capture known biology using this balancing approach compared to just applying Pearson correlation. They look to see if groups of related genes are from the same GO term (i.e. do genes from the same GO term cluster together?)
#
# 2. From [Zhu et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4768301/), we apply their "gene hubbiness correction" procedure. This is part of their [SEEK](https://seek.princeton.edu/seek/) algorithm. This correction procedure is motivated by the observation that hubby or well-connected genes in the co-expression network represent global, well-co-expressed processes, and can contaminate the search results regardless of query composition due to the effect of unbalanced gene connectivity in a scale-free co-expression network, and can lead to non-specific results in search or clustering approaches. To avoid the bias created by hubby genes they correct for each gene g’s correlation by subtracting g's average correlation.

# %load_ext autoreload
# %autoreload 2
import os
import pandas as pd
import plotnine as pn
import seaborn as sns
import matplotlib.pyplot as plt
import umap
import random
import numpy as np
from sklearn.decomposition import PCA
from core_acc_modules import paths

# ## Set user parameters
#
# For now we will vary the correlation threshold (`corr_threshold`) but keep the other parameters consistent
#
# We will run this notebook for each threshold parameter

# Params
num_PCs = 20

# Load expression data
pao1_compendium_filename = paths.PAO1_COMPENDIUM
pa14_compendium_filename = paths.PA14_COMPENDIUM

pao1_compendium = pd.read_csv(pao1_compendium_filename, sep="\t", header=0, index_col=0)
pa14_compendium = pd.read_csv(pa14_compendium_filename, sep="\t", header=0, index_col=0)

print(pao1_compendium.shape)
pao1_compendium.head()

print(pa14_compendium.shape)
pa14_compendium.head()

# ## Correlation of raw gene expression data
#
# Use this as a comparison to see how the correlations are changed after each correction method

# Correlation
pao1_corr_original = pao1_compendium.corr()
pa14_corr_original = pa14_compendium.corr()

pao1_corr_original.head(10)

pa14_corr_original.head(10)

# ## SPELL + correlation
#
# _Review of SVD_
#
# Singular Value Decomposition is a way to factorize your matrix, $X^{mxn}$ into singular vectors and singular values: $X = U \Sigma V^*$
#
# In our case $X$ is **gene x sample** and then the columns of $U$ (gene x gene) are the left singular vectors (gene coefficient vectors); $\Sigma$ (the same dimensions as $X$) has singular values and is diagonal (mode amplitudes); and $V^T$ (sample x sample) has rows that are the right singular vectors (expression level vectors).
#
# Here we are using SVD to reduce the noise in our original data by performing dimensionality reduction. This dimensionality is done by neglecting the small singular values in the diagonal matrix $\Sigma$. Normally people would get the projection of the original data onto the singular vectors by $U \Sigma$ and apply the correlation on the projected data. Here, we're following the description in [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926) where they performed correlation on $U$ only.

# Transpose compendia to be gene x sample
# Here we're interested in how genes cluster
pao1_compendium_T = pao1_compendium.T
pa14_compendium_T = pa14_compendium.T

# Apply SVD
pao1_U, pao1_s, pao1_Vh = np.linalg.svd(pao1_compendium_T, full_matrices=False)
pa14_U, pa14_s, pa14_Vh = np.linalg.svd(pa14_compendium_T, full_matrices=False)

print(pao1_compendium.shape)
print(pao1_U.shape, pao1_s.shape, pao1_Vh.shape)

# In the graph, we can see that although we have 847 singular values in s, most of those (after the 20th entry or so) are pretty small. So it might make sense to use only the information related to the first (say, 20) singular values to build a reduced representation.

plt.plot(pao1_s[:100])

print(pa14_compendium.shape)
print(pa14_U.shape, pa14_s.shape, pa14_Vh.shape)

# Convert ndarray to df to use corr()
pao1_U_df = pd.DataFrame(
    data=pao1_U,
    index=pao1_compendium_T.index,  # columns=pao1_compendium_T.index
)
pa14_U_df = pd.DataFrame(
    data=pa14_U,
    index=pa14_compendium_T.index,  # columns=pa14_compendium_T.index
)

pao1_U_df.head()

# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
num_singular_values = 20
pao1_corr_spell = pao1_U_df.iloc[:, :num_singular_values].T.corr()
pa14_corr_spell = pa14_U_df.iloc[:, :num_singular_values].T.corr()

# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_corr_spell.abs(), cmap="viridis")
h1.fig.suptitle(
    f"Correlation of PAO1 genes (SPELL corrected using {num_singular_values} vectors)",
    y=1.05,
)

# Would expect SVD and PCA to dappen highlight high correlation signals, but this is not what we see. At least this issue seems to be consistent between SVD and PCA. So either
# * Implemented both incorrectly
# * This is real but what explains this pattern

plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_corr_spell.abs(), cmap="viridis")
h2.fig.suptitle(
    f"Correlation of PA14 genes (SPELL corrected using {num_singular_values} vectors)",
    y=1.05,
)

# ## PCA + correlation

# +
# Embed expression data into low dimensional space
pca = PCA(n_components=num_PCs)

model_pao1 = pca.fit(pao1_compendium_T)
pao1_encoded = model_pao1.transform(pao1_compendium_T)

pao1_encoded_df = pd.DataFrame(data=pao1_encoded, index=pao1_compendium_T.index)

print(pao1_encoded_df.shape)

# +
model_pa14 = pca.fit(pa14_compendium_T)
pa14_encoded = model_pa14.transform(pa14_compendium_T)

pa14_encoded_df = pd.DataFrame(data=pa14_encoded, index=pa14_compendium_T.index)
print(pa14_encoded_df.shape)
# -

# Correlation
pao1_corr_pca = pao1_encoded_df.T.corr()
pa14_corr_pca = pa14_encoded_df.T.corr()

# Plot heatmap
plt.figure(figsize=(20, 20))
h3 = sns.clustermap(pao1_corr_pca.abs(), cmap="viridis")
h3.fig.suptitle(f"Correlation of {num_PCs} PCA encoded PAO1 genes", y=1.05)

plt.figure(figsize=(20, 20))
h4 = sns.clustermap(pa14_corr_pca.abs(), cmap="viridis")
h4.fig.suptitle(f"Correlation of {num_PCs} PCA encoded PA14 genes", y=1.05)

# ## SEEK
#
# Look into https://het.io/ correction

"""# Correlation of U
pao1_corr_seek = pao1_compendium.corr()
pa14_corr_seek = pa14_compendium.corr()"""

"""# Subtract mean correlation score for gene
pao1_corr_mean = pao1_corr_seek.mean()
pa14_corr_mean = pa14_corr_seek.mean()

pao1_corr_mean"""

"""pao1_corr_corrected = pao1_corr_seek - pao1_corr_mean
pa14_corr_corrected = pa14_corr_seek - pa14_corr_mean"""

# +
# pao1_corr_seek.head()

# +
# Notice that the corrected correlation matrix is not symmetric anymore because we
# are subtracting the mean for each row
# pao1_corr_corrected.head()

# +
# pa14_corr_corrected.head()

# +
# Create a similarity matrix usingn the threshold defined above
# The similarity matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by by the correlation score, otherwise the value is set to 0
# pao1_corr_corrected[pao1_corr_corrected.abs() < corr_threshold] = 0.0
# pa14_corr_corrected[pa14_corr_corrected.abs() < corr_threshold] = 0.0

# pao1_corr_corrected.head()
# -

"""# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_corr_corrected.abs(), cmap="viridis")
h1.fig.suptitle(f"Correlation of PAO1 genes using (SEEK corrected)")

# Save
pao1_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_corr_SEEK_clustermap.png"
)
# h1.savefig(pao1_clustermap_filename, dpi=300)"""

"""plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_corr_corrected.abs(), cmap="viridis")
h2.fig.suptitle(f"Correlation of PA14 genes (SEEK corrected)")

# Save
pa14_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_corr_SEEK_clustermap.png"
)
# h2.savefig(pa14_clustermap_filename, dpi=300)"""

# ## Hetio

# **Takeaway:**
