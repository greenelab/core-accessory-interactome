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
# When we performed clustering on the correlation matrices (using Pearson correlation) we found that pairs of genes had either very high correlation scores (>0.5) or very low correlation scores (<0.1). As a result gene pairs that were highly correlated clustered into a single large module. This finding is consistent with a [previous study](https://link.springer.com/article/10.1186/1471-2164-7-187), which found that KEGG (a database that containes genes or proteins annotated with specific biological processes as reported in the literature) is bias in some biological processes represented. Figure 1C demonstrates that a large fraction of gene pairs are ribosomal relationships - in the top 0.1% most co-expressed genes, 99% belong to the ribosome pathway.
# Furthermore, protein function prediction based on co-expression drop dramatically after removing the ribisome pathway (Figure 1A, B).
#
# This notebook applies different corrections to try to remove this very dominant global signal in the data. This notebook follows from [1a_transformation_correlation_analysis.ipynb](1a_transformation_correlation_analysis.ipynb). Here we are applying dimensionality reduction techniques in addition to scaling the data.

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
# Here we set the number of PCs or singular vectors to use. We are starting with 300 since this is what [eADAGE](https://pubmed.ncbi.nlm.nih.gov/28711280/) used.

# Params
num_PCs = 300
num_PCs_log = 100
num_singular_values = 300
num_singular_values_log = 100

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
# Here is the correlation of the raw data without any malnipulations. This will serve as a reference to compare the correlations below where applied corrections to the correlations to account for the dominant signal described above.

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
# Here I am trying to create a null distribution of correlation score to determine if a pair of genes are more correlated than expected. This will still require some thought, but for now this is a first pass.

# Shuffle values per gene
pao1_shuffled_compendium = pao1_compendium.apply(lambda x: x.sample(frac=1).values)
pa14_shuffled_compendium = pa14_compendium.apply(lambda x: x.sample(frac=1).values)

"""# Correlation
pao1_corr_shuffled = pao1_shuffled_compendium.corr()
pa14_corr_shuffled = pa14_shuffled_compendium.corr()"""

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
# From [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926), we apply their "signal balancing technique that enhances biological information". This is the first part of their [SPELL](https://spell.yeastgenome.org/) algorithm that is described in section 2.3.1. SPELL calculates the correlation on the gene coefficient matrix, $U$ (i.e. how much genes contribute to a latent variable) that is generated after applying SVD. This matrix represents how genes contribute to independent latent variables that capture the signal in the data where the variance of the variables is 1. The idea is that correlations between gene contributions are more balanced so that less prominent patterns are amplified and more dominant patterns are dampended due to this compression. Figure 3 shows how well SPELL recapitulates biology (i.e. the relationship between genes within a GO term) compared to Pearson correlation.

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

# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
pao1_corr_spell = pao1_U_df.iloc[:, :num_singular_values].T.corr()
pa14_corr_spell = pa14_U_df.iloc[:, :num_singular_values].T.corr()

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
pao1_compendium_log10 = np.log10(1 + pao1_compendium_T)
pa14_compendium_log10 = np.log10(1 + pa14_compendium_T)

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

"""# Plot heatmap
plt.figure(figsize=(20, 20))
h1a = sns.clustermap(pao1_corr_spell.abs(), cmap="viridis")
h1a.fig.suptitle(
    f"log transform + SPELL corrected using {num_singular_values} vectors (PAO1)",
    y=1.05,
)"""

"""plt.figure(figsize=(20, 20))
h2a = sns.clustermap(pa14_corr_spell.abs(), cmap="viridis")
h2a.fig.suptitle(
    f"log transformed + SPELL corrected using {num_singular_values} vectors (PA14)",
    y=1.05,
)"""

# ## PCA + correlation
#
# Here we are going to calculate the correlation of the gene coefficients, stored in the PC components matrix. We expect this correlation matrix to look very similar to the SPELL results. This is acting as a positive control that we implemented the SPELL correction correctly

# +
# Embed expression data into low dimensional space
pca = PCA(n_components=num_PCs)

model_pao1 = pca.fit(pao1_compendium)

# Using components matrix: gene x PCs (each PC value is how much a gene contributes)
print(model_pao1.components_.shape)
pao1_pc_weights_df = pd.DataFrame(
    data=model_pao1.components_, columns=pao1_compendium.columns
)
pao1_pc_weights_df.head()

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

# Using components matrix: gene x PCs (each PC value is how much a gene contributes)
print(model_pa14.components_.shape)
pa14_pc_weights_df = pd.DataFrame(
    data=model_pa14.components_, columns=pa14_compendium.columns
)

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
pao1_corr_pca = pao1_pc_weights_df.corr()
pa14_corr_pca = pa14_pc_weights_df.corr()

"""# Plot heatmap
plt.figure(figsize=(20, 20))
h3 = sns.clustermap(pao1_corr_pca.abs(), cmap="viridis")
h3.fig.suptitle(f"PCA using {num_PCs} PCs (PAO1 genes)", y=1.05)"""

"""plt.figure(figsize=(20, 20))
h4 = sns.clustermap(pa14_corr_pca.abs(), cmap="viridis")
h4.fig.suptitle(f"PCA using {num_PCs} PCs (PA14 genes)", y=1.05)"""

# ## log transform + PCA + correlation

# log transform data
# Note: add 1 to avoid -inf and so 0 maps to those
pao1_compendium_log10 = np.log10(1 + pao1_compendium)
pa14_compendium_log10 = np.log10(1 + pa14_compendium)

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

"""# Plot heatmap
plt.figure(figsize=(20, 20))
h3a = sns.clustermap(pao1_corr_pca.abs(), cmap="viridis")
h3a.fig.suptitle(f"log transformed + PCA using {num_PCs} PCs (PAO1 genes)", y=1.05)"""

"""plt.figure(figsize=(20, 20))
h4a = sns.clustermap(pa14_corr_pca.abs(), cmap="viridis")
h4a.fig.suptitle(f"log transformed + PCA using {num_PCs} PCs (PA14 genes)", y=1.05)"""

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

"""# Plot heatmap
plt.figure(figsize=(20, 20))
h5 = sns.clustermap(pao1_corr_seek.abs(), cmap="viridis")
h5.fig.suptitle("Correlation of PAO1 genes using (SEEK corrected)", y=1.05)"""

"""plt.figure(figsize=(20, 20))
h6 = sns.clustermap(pa14_corr_seek.abs(), cmap="viridis")
h6.fig.suptitle("Correlation of PA14 genes (SEEK corrected)", y=1.05)"""

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

"""# Plot heatmap
plt.figure(figsize=(20, 20))
h7 = sns.clustermap(pao1_corr_hetio_df.abs(), cmap="viridis")
h7.fig.suptitle("Correlation of PAO1 genes using (Hetio corrected)", y=1.05)"""

"""# Plot heatmap
plt.figure(figsize=(20, 20))
h8 = sns.clustermap(pa14_corr_hetio_df.abs(), cmap="viridis")
h8.fig.suptitle("Correlation of PA14 genes using (Hetio corrected)", y=1.05)"""

# ## Subtract the mean
#
# Inspired by hetio search, here we are subtracting the mean of the row and column from each correlation value.
#
# _implementation note:_
#
# We are subtracting the row mean and column mean separately based on a simple proof where $A_ij$ is the value in row _i_ and column _j_, There are a total of _N_ columns and rows since the matrix is square.
#
# $A_{ij} - \frac{\Sigma_{x=0}^N A_{xj} - \Sigma_{x=0}^N A_{ix}}{2N}$
#
# $A_{ij} - \frac{\Sigma_{x=0}^N A_{xj}}{2N} - \frac{\Sigma_{x=0}^N A_{ix}}{2N}$
#
# $A_{ij} - \frac{1}{2}\frac{\Sigma_{x=0}^N A_{xj}}{N} - \frac{1}{2}\frac{\Sigma_{x=0}^N A_{ix}}{N}$

pao1_corr_mean_subtract = pao1_corr_original.subtract(
    pao1_corr_original.mean(axis=0) / 2, axis=0
).subtract(pao1_corr_original.mean(axis=1) / 2, axis=1)

pa14_corr_mean_subtract = pa14_corr_original.subtract(
    pa14_corr_original.mean(axis=0) / 2, axis=0
).subtract(pa14_corr_original.mean(axis=1) / 2, axis=1)

pao1_corr_mean_subtract.head()

# Plot heatmap
plt.figure(figsize=(20, 20))
h9 = sns.clustermap(pao1_corr_mean_subtract.abs(), cmap="viridis")
h9.fig.suptitle("Subtract mean from correlation of PAO1 genes", y=1.05)

# Plot heatmap
plt.figure(figsize=(20, 20))
h10 = sns.clustermap(pa14_corr_mean_subtract.abs(), cmap="viridis")
h10.fig.suptitle("Subtract mean from correlation of PA14 genes", y=1.05)

# **Takeaway:**
#
# * PCA vs SVD
#
#     Given $ X W = Z$.
#
#     The goal of **Principal Component Analysis (PCA)** is to find a weight matrix ($W$) that reduces the data into a low dimensional space that captures the as much information from our input data, $X$, as possible. In other words, we want to find a $W$ that captures as much of the variance in the original data as possible.
#
#     We can use the covariance matrix to describe the input data $X$. $Cov(X)$ is a symmetric matrix that has variances along the diagonal (i.e. spread of the data) and covariances on the off diagonal (orientation of the data).
#
#     We can factorize $Cov(X) = VDV^T$, where $V$ are the eigenvectors, which represent the direction of variance in our data. $V$ are the our principal components to project our data onto. So the weight matrix is composed of these principal components. And when you multiple $XW = XVD = Z$, which is our data projected onto the most variable directions
#
#     **Singular Value Decomposition (SVD)** is a way to factorize your matrix, $X^{mxn}$ into singular vectors and singular values: $X = U \Sigma V^*$
#
#     In our case $X$ is **gene x sample** and then the columns of $U$ are the left singular vectors; $\Sigma$ (eigengene x eigensample) has singular values and is diagonal (mode amplitudes); and $V^T$ has rows that are the right singular vectors.
#
#     Overall, PCA is performing SVD but using the covariance matrix as input. Using the covariance matrix as the input means that the data is centered (i.e. mean is 0). So there is no difference and we might want to use SVD in this case because there is an associated publication to cite.
#
#
# * log transform will compress the data so that highly variable genes don't dominant the correlation.
#
# * Subtracting the mean correlation score is trying to account for a baseline signal that is present in the data. This still results in fairly large clusters.
#
# Based on this exploration, we will use log transform + SVD correlation matrix to identify clusters.
