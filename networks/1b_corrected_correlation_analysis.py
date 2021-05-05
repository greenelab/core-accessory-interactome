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
from sklearn import preprocessing
import matplotlib.pyplot as plt
import umap
import random
import numpy as np
import scipy
from sklearn.decomposition import PCA
from core_acc_modules import paths

# ## Set user parameters
#
# For now we will vary the correlation threshold (`corr_threshold`) but keep the other parameters consistent
#
# We will run this notebook for each threshold parameter

# Params
num_PCs = 100
num_singular_values = 100
corr_threshold = 0.5

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

pao1_corr_original.head()

pa14_corr_original.head()

"""# Plot heatmap
plt.figure(figsize=(20, 20))
o1 = sns.clustermap(pao1_corr_original.abs(), cmap="viridis")
o1.fig.suptitle(f"Correlation of raw PAO1 genes", y=1.05)"""

"""# Plot heatmap
plt.figure(figsize=(20, 20))
o2 = sns.clustermap(pa14_corr_original.abs(), cmap="viridis")
o2.fig.suptitle(f"Correlation of raw PA14 genes", y=1.05)"""

# ## Correlation of permuted raw gene expression data
#

# Shuffle values per gene
pao1_shuffled_compendium = pao1_compendium.apply(lambda x: x.sample(frac=1).values)
pa14_shuffled_compendium = pa14_compendium.apply(lambda x: x.sample(frac=1).values)

# Correlation
pao1_corr_shuffled = pao1_shuffled_compendium.corr()
pa14_corr_shuffled = pa14_shuffled_compendium.corr()

"""# Plot heatmap
plt.figure(figsize=(20, 20))
o3 = sns.clustermap(pao1_corr_shuffled.abs(), cmap="viridis")
o3.fig.suptitle(f"Correlation of Shuffled PAO1 genes", y=1.05)"""

"""# Plot heatmap
plt.figure(figsize=(20, 20))
o4 = sns.clustermap(pa14_corr_shuffled.abs(), cmap="viridis")
o4.fig.suptitle(f"Correlation of Shuffled PA14 genes", y=1.05)"""

# ## SPELL + correlation
#
# _Review of SVD_
#
# Singular Value Decomposition is a way to factorize your matrix, $X^{mxn}$ into singular vectors and singular values: $X = U \Sigma V^*$
#
# In our case $X$ is **gene x sample** and then the columns of $U$ (gene x eigensample) are the left singular vectors (gene coefficient vectors); $\Sigma$ (eigengene x eigensample) has singular values and is diagonal (mode amplitudes); and $V^T$ (eigengene x sample) has rows that are the right singular vectors (expression level vectors).
#
# Here we are using SVD to reduce the noise in our original data by performing dimensionality reduction. This dimensionality is done by neglecting the small singular values in the diagonal matrix $\Sigma$. Normally people would get the projection of the original data onto the singular vectors by $U \Sigma$ and apply the correlation on the projected data. Here, we're following the description in [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926) where they performed correlation on $U$ only.
#
# From [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926), we apply their "signal balancing technique that enhances biological information". This is the first part of their [SPELL](https://spell.yeastgenome.org/) algorithm that is described in section 2.3.1. They apply SVD to the gene expression matrix to reduce the noise that can lead to spurious results and then apply correlation between genes using the singular vector.  Correlations between genes in U equally weight each dimension of the orthonormal basis and balance their contributions such that the least prominent patterns are amplified and more dominant patterns are dampened. This process helps reveal biological signals, as some of the dominant patterns in many microarray datasets are not biologically meaningful. These dominant signals might be due to technical artifacts that lead to artificially high correlation levels (i.e. If a set of samples have all genes highly expressed or lowly expressed). Section 2.4 evaluates SPELL on its ability to capture known biology using this balancing approach compared to just applying Pearson correlation. They look to see if groups of related genes are from the same GO term (i.e. do genes from the same GO term cluster together?)

# Transpose compendia to be gene x sample
# Here we're interested in how genes cluster
pao1_compendium_T = pao1_compendium.T
pa14_compendium_T = pa14_compendium.T

# Apply SVD
pao1_U, pao1_s, pao1_Vh = np.linalg.svd(pao1_compendium_T, full_matrices=False)
pa14_U, pa14_s, pa14_Vh = np.linalg.svd(pa14_compendium_T, full_matrices=False)

print(pao1_compendium_T.shape)
print(pao1_U.shape, pao1_s.shape, pao1_Vh.shape)

# In the graph, we can see that although we have 847 singular values in s, most of those (after the 20th entry or so) are pretty small. So it might make sense to use only the information related to the first (say, 20) singular values to build a reduced representation.

plt.plot(pao1_s[:100])

print(pa14_compendium_T.shape)
print(pa14_U.shape, pa14_s.shape, pa14_Vh.shape)

# Convert ndarray to df to use corr()
pao1_U_df = pd.DataFrame(
    data=pao1_U,
    index=pao1_compendium_T.index,  # columns=pao1_compendium.index
)
pa14_U_df = pd.DataFrame(
    data=pa14_U,
    index=pa14_compendium_T.index,  # columns=pa14_compendium.index
)

pao1_U_df.head()

# +
# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
pao1_corr_spell = pao1_U_df.iloc[:, :num_singular_values].T.corr()
pa14_corr_spell = pa14_U_df.iloc[:, :num_singular_values].T.corr()

# pao1_corr_spell = pao1_U_df.T.corr()
# pa14_corr_spell = pa14_U_df.T.corr()
# -

# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_corr_spell.abs(), cmap="viridis")
h1.fig.suptitle(
    f"Correlation of PAO1 genes (SPELL corrected using {num_singular_values} vectors)",
    y=1.05,
)

plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_corr_spell.abs(), cmap="viridis")
h2.fig.suptitle(
    f"Correlation of PA14 genes (SPELL corrected using {num_singular_values} vectors)",
    y=1.05,
)

# ## log transform + SPELL + correlation

# log transform data
pao1_compendium_log10 = np.log10(pao1_compendium_T)
pa14_compendium_log10 = np.log10(pa14_compendium_T)

# Set inf to 0
pao1_compendium_log10[np.isinf(pao1_compendium_log10)] = 0
pa14_compendium_log10[np.isinf(pa14_compendium_log10)] = 0

# Apply SVD
pao1_U, pao1_s, pao1_Vh = np.linalg.svd(pao1_compendium_log10, full_matrices=False)
pa14_U, pa14_s, pa14_Vh = np.linalg.svd(pa14_compendium_log10, full_matrices=False)

print(pao1_compendium_T.shape)
print(pao1_U.shape, pao1_s.shape, pao1_Vh.shape)

print(pa14_compendium_T.shape)
print(pa14_U.shape, pa14_s.shape, pa14_Vh.shape)

# Convert ndarray to df to use corr()
pao1_U_df = pd.DataFrame(
    data=pao1_U,
    index=pao1_compendium_T.index,  # columns=pao1_compendium.index
)
pa14_U_df = pd.DataFrame(
    data=pa14_U,
    index=pa14_compendium_T.index,  # columns=pa14_compendium.index
)

# +
# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U

pao1_corr_spell = pao1_U_df.iloc[:, :num_singular_values].T.corr()
pa14_corr_spell = pa14_U_df.iloc[:, :num_singular_values].T.corr()
# -

# Plot heatmap
plt.figure(figsize=(20, 20))
h1a = sns.clustermap(pao1_corr_spell.abs(), cmap="viridis")
h1a.fig.suptitle(
    f"log transform + SPELL corrected using {num_singular_values} vectors (PAO1)",
    y=1.05,
)

plt.figure(figsize=(20, 20))
h2a = sns.clustermap(pa14_corr_spell.abs(), cmap="viridis")
h2a.fig.suptitle(
    f"log transformed + SPELL corrected using {num_singular_values} vectors (PA14)",
    y=1.05,
)

# ## PCA + correlation
#
# Here we try applying the correlation to the data matrix projected onto the reduced space as well as looking at the correlation on the components matrix that contains gene coefficients per PC.

# +
# Embed expression data into low dimensional space
pca = PCA(n_components=num_PCs)

model_pao1 = pca.fit(pao1_compendium)
pao1_encoded = model_pao1.transform(pao1_compendium)

# Using components matrix: gene x PCs (each PC value is how much a gene contributes)
print(model_pao1.components_.shape)
pao1_pc_weights_df = pd.DataFrame(
    data=model_pao1.components_, columns=pao1_compendium.columns
)
pao1_pc_weights_df.head()

# Using reduced matrix: gene x PCs (each PC is a linear combination of samples)
# pao1_encoded_df = pd.DataFrame(data=pao1_encoded, index=pao1_compendium_T.index)
# print(pao1_encoded_df.shape)

# +
pca_variance = model_pao1.explained_variance_

plt.figure(figsize=(8, 6))
plt.bar(
    range(num_PCs), pca_variance, alpha=0.5, align="center", label="individual variance"
)
plt.legend()
plt.title("Variance explained for PAO1")
plt.ylabel("Variance ratio")
plt.xlabel("Principal components")
plt.show()

# +
model_pa14 = pca.fit(pa14_compendium)
pa14_encoded = model_pa14.transform(pa14_compendium)

# Using components matrix: gene x PCs (each PC value is how much a gene contributes)
print(model_pa14.components_.shape)
pa14_pc_weights_df = pd.DataFrame(
    data=model_pa14.components_, columns=pa14_compendium.columns
)

# Using reduced matrix: gene x PCs (each PC is a linear combination of samples)
# pa14_encoded_df = pd.DataFrame(data=pa14_encoded, index=pa14_compendium.index)
# print(pa14_encoded_df.shape)

# +
pca_variance = model_pa14.explained_variance_

plt.figure(figsize=(8, 6))
plt.bar(
    range(num_PCs), pca_variance, alpha=0.5, align="center", label="individual variance"
)
plt.legend()
plt.title("Variance explained for PA14")
plt.ylabel("Variance ratio")
plt.xlabel("Principal components")
plt.show()
# -

# Correlation
# pao1_corr_pca = pao1_encoded_df.T.corr()
# pa14_corr_pca = pa14_encoded_df.T.corr()
pao1_corr_pca = pao1_pc_weights_df.corr()
pa14_corr_pca = pa14_pc_weights_df.corr()

# Plot heatmap
plt.figure(figsize=(20, 20))
h3 = sns.clustermap(pao1_corr_pca.abs(), cmap="viridis")
h3.fig.suptitle(f"PCA using {num_PCs} PCs (PAO1 genes)", y=1.05)

plt.figure(figsize=(20, 20))
h4 = sns.clustermap(pa14_corr_pca.abs(), cmap="viridis")
h4.fig.suptitle(f"PCA using {num_PCs} PCs (PA14 genes)", y=1.05)

# Would expect SVD and PCA to dappen highlight high correlation signals, but this is not what we see. At least this issue seems to be consistent between SVD and PCA. So either
# * Implemented both incorrectly
# * This is real but what explains this pattern

# ## log transform + PCA + correlation

# log transform data
pao1_compendium_log10 = np.log10(pao1_compendium)
pa14_compendium_log10 = np.log10(pa14_compendium)

# Set inf to 0
pao1_compendium_log10[np.isinf(pao1_compendium_log10)] = 0
pa14_compendium_log10[np.isinf(pa14_compendium_log10)] = 0

# +
# Embed expression data into low dimensional space
pca = PCA(n_components=num_PCs)

model_pao1 = pca.fit(pao1_compendium_log10)

# Using components matrix: gene x PCs (each PC value is how much a gene contributes)
print(model_pao1.components_.shape)
pao1_pc_weights_df = pd.DataFrame(
    data=model_pao1.components_, columns=pao1_compendium_log10.columns
)
pao1_pc_weights_df.head()


# +
pca_variance = model_pao1.explained_variance_

plt.figure(figsize=(8, 6))
plt.bar(
    range(num_PCs), pca_variance, alpha=0.5, align="center", label="individual variance"
)
plt.legend()
plt.title("Variance explained for log transformed PAO1")
plt.ylabel("Variance ratio")
plt.xlabel("Principal components")
plt.show()

# +
model_pa14 = pca.fit(pa14_compendium_log10)
pa14_encoded = model_pa14.transform(pa14_compendium_log10)

# Using components matrix: gene x PCs (each PC value is how much a gene contributes)
print(model_pa14.components_.shape)
pa14_pc_weights_df = pd.DataFrame(
    data=model_pa14.components_, columns=pa14_compendium_log10.columns
)

# +
pca_variance = model_pa14.explained_variance_

plt.figure(figsize=(8, 6))
plt.bar(
    range(num_PCs), pca_variance, alpha=0.5, align="center", label="individual variance"
)
plt.legend()
plt.title("Variance explained for log transformed PA14")
plt.ylabel("Variance ratio")
plt.xlabel("Principal components")
plt.show()
# -

# Correlation
pao1_corr_pca = pao1_pc_weights_df.corr()
pa14_corr_pca = pa14_pc_weights_df.corr()

# Plot heatmap
plt.figure(figsize=(20, 20))
h3a = sns.clustermap(pao1_corr_pca.abs(), cmap="viridis")
h3a.fig.suptitle(f"log transformed + PCA using {num_PCs} PCs (PAO1 genes)", y=1.05)

plt.figure(figsize=(20, 20))
h4a = sns.clustermap(pa14_corr_pca.abs(), cmap="viridis")
h4a.fig.suptitle(f"log transformed + PCA using {num_PCs} PCs (PA14 genes)", y=1.05)

# ## SEEK + correlation
#
# From [Zhu et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4768301/), we apply their "gene hubbiness correction" procedure. This is part of their [SEEK](https://seek.princeton.edu/seek/) algorithm. This correction procedure is motivated by the observation that hubby or well-connected genes in the co-expression network represent global, well-co-expressed processes, and can contaminate the search results regardless of query composition due to the effect of unbalanced gene connectivity in a scale-free co-expression network, and can lead to non-specific results in search or clustering approaches. To avoid the bias created by hubby genes they correct for each gene g’s correlation by subtracting g's average correlation.
#
# We noticed that the corrected correlation matrix is not symmetric anymore because we are subtracting the mean for each row. Instead we will normalize the correlation scores per gene

# +
# 0-1 normalize per gene
scaler = preprocessing.MinMaxScaler()

pao1_corr_scaled = scaler.fit_transform(pao1_corr_original)
pa14_corr_scaled = scaler.fit_transform(pa14_corr_original)

pao1_corr_seek = pd.DataFrame(
    pao1_corr_scaled, columns=pao1_corr_original.columns, index=pao1_corr_original.index
)
pa14_corr_seek = pd.DataFrame(
    pa14_corr_scaled, columns=pa14_corr_original.columns, index=pa14_corr_original.index
)
# -

pao1_corr_original.head()

pao1_corr_seek.head()

# Plot heatmap
plt.figure(figsize=(20, 20))
h5 = sns.clustermap(pao1_corr_seek.abs(), cmap="viridis")
h5.fig.suptitle("Correlation of PAO1 genes using (SEEK corrected)")

plt.figure(figsize=(20, 20))
h6 = sns.clustermap(pa14_corr_seek.abs(), cmap="viridis")
h6.fig.suptitle("Correlation of PA14 genes (SEEK corrected)")

# ## Hetio
#
# https://het.io/ allows users to search a heterogeneous network (i.e. are graphs with multiple node and edge types like disease or gene nodes) to find connections between entities such as eye and breast cancer. It is used as a hypothesis generating tool.
#
# To search the network, they empoloy a correction for hubbiness. Here we use their correction. This correction is similar to SEEK in that it is applying a correction using both the source and target node to deal with the asymmetry that we encountered.
#
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4497619/ see "Feature computation metrics" in Methods section:
# we developed the degree-weighted path count (DWPC) which individually downweights each path between a source and target node. Each path receives a path-degree product (PDP) calculated by: 1) extracting all metaedge-specific degrees along the path (D path), where each edge composing the path contributes two degrees; 2) raising each degree to the −w power, where w ≥ 0 and is called the damping exponent; 3) multiplying all exponentiated degrees to yield the PDP.

# Start with an adjacency matrix
pao1_adj = (pao1_corr_original.abs() >= corr_threshold).astype(float)
pa14_adj = (pa14_corr_original.abs() >= corr_threshold).astype(float)

pao1_matrix = np.array(pao1_adj)
pa14_matrix = np.array(pa14_adj)

# Sum along the columns to get the degree per node
pao1_vector = np.array(pao1_matrix.sum(axis=0)).flatten()
pa14_vector = np.array(pa14_matrix.sum(axis=0)).flatten()

# +
# Multiply by scaled sum vector
damping_exponent = 0.1
with np.errstate(divide="ignore"):
    # Scale node degree
    pao1_vector **= -damping_exponent
    pa14_vector **= -damping_exponent

# Set inf to 0
pao1_vector[np.isinf(pao1_vector)] = 0
pa14_vector[np.isinf(pa14_vector)] = 0

# Create a matrix with damping factor on the diagonals
pao1_vector = scipy.sparse.diags(pao1_vector)
pa14_vector = scipy.sparse.diags(pa14_vector)
# -

# matrix multiplication
pao1_corr_hetio = pao1_matrix @ pao1_vector
pa14_corr_hetio = pa14_matrix @ pa14_vector

pao1_corr_hetio_df = pd.DataFrame(data=pao1_corr_hetio)
pa14_corr_hetio_df = pd.DataFrame(data=pa14_corr_hetio)

pao1_corr_hetio_df.head()

pa14_corr_hetio_df.head()

# Plot heatmap
plt.figure(figsize=(20, 20))
h7 = sns.clustermap(pao1_corr_hetio_df.abs(), cmap="viridis")
h7.fig.suptitle("Correlation of PAO1 genes using (Hetio corrected)")

# Plot heatmap
plt.figure(figsize=(20, 20))
h8 = sns.clustermap(pa14_corr_hetio_df.abs(), cmap="viridis")
h8.fig.suptitle("Correlation of PA14 genes using (Hetio corrected)")

# ## Subtract the mean

# **Takeaway:**
