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
# 1. From [Hibbs et. al.](https://academic.oup.com/bioinformatics/article/23/20/2692/229926), we apply their "signal balancing technique that enhances biological information". This is the first part of their [SPELL](https://spell.yeastgenome.org/) algorithm that is described in section 2.3.1. They apply SVD to the gene expression matrix to reduce the noise that can lead to spurious results and then apply correlation between genes using the singular vector.  Correlations between genes in U equally weight each dimension of the orthonormal basis and balance their contributions such that the least prominent patterns are amplified and more dominant patterns are dampened. This process helps reveal biological signals, as some of the dominant patterns in many microarray datasets are not biologically meaningful (**WHAT ARE THESE DOMINANT SIGNALS???**) Section 2.4 evaluates SPELL on its ability to capture known biology using this balancing approach compared to just applying Pearson correlation. They look to see if groups of related genes are from the same GO term (i.e. do genes from the same GO term cluster together?)
#
# 2. From [Zhu et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4768301/), we apply their "gene hubbiness correction" procedure. This is part of their [SEEK] algorithm. This correction procedure is motivated by the observation that hubby or well-connected genes in the co-expression network represent global, well-co-expressed processes, and can contaminate the search results regardless of query composition due to the effect of unbalanced gene connectivity in a scale-free co-expression network, and can lead to non-specific results in search or clustering approaches. To avoid the bias created by hubby genes they correct for each gene g’s correlation by subtracting g's average correlation.

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
from scipy.spatial.distance import pdist, squareform
from core_acc_modules import paths

# ## Set user parameters
#
# For now we will vary the correlation threshold (`corr_threshold`) but keep the other parameters consistent
#
# We will run this notebook for each threshold parameter

# +
# Params
corr_threshold = 0.5

# Output files
pao1_membership_filename = f"pao1_membership_{corr_threshold}.tsv"
pa14_membership_filename = f"pa14_membership_{corr_threshold}.tsv"
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

# ## Correlation of raw gene expression data
#
# Use this as a comparison to see how the correlations are changed after each correction method

# Correlation
pao1_corr_original = pao1_compendium.corr()
pa14_corr_original = pa14_compendium.corr()

pao1_corr_original.head(10)

pa14_corr_original.head(10)

# ## SPELL
#
# _Review of SVD_
#
# Singular Value Decomposition is a way to decompose your matrix, $X^{mxn}$ into 3 matrices: $X = U \Sigma V^*$
#
# In our case $X$ is gene x sample and then the columns of $U$ (gene x gene) are the left singular vectors (gene coefficient vectors); $\Sigma$ (the same dimensions as $X$) has singular values and is diagonal (mode amplitudes); and $V^T$ (sample x sample) has rows that are the right singular vectors (expression level vectors).

# Transpose compendia to be gene x sample
pao1_compendium_T = pao1_compendium.T
pa14_compendium_T = pa14_compendium.T

# Apply SVD
pao1_U, pao1_s, pao1_Vh = np.linalg.svd(pao1_compendium_T, full_matrices=False)
pa14_U, pa14_s, pa14_Vh = np.linalg.svd(pa14_compendium_T, full_matrices=False)

print(pao1_compendium.shape)
print(pao1_U.shape, pao1_s.shape, pao1_Vh.shape)

print(pa14_compendium.shape)
print(pa14_U.shape, pa14_s.shape, pa14_Vh.shape)

# Convert ndarray to df to use corr()
pao1_U_df = pd.DataFrame(
    data=pao1_U, index=pao1_compendium_T.index, columns=pao1_compendium_T.columns
)
pa14_U_df = pd.DataFrame(
    data=pa14_U, index=pa14_compendium_T.index, columns=pa14_compendium_T.columns
)

pao1_U_df.head()

# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
pao1_corr_spell = pao1_U_df.T.corr()
pa14_corr_spell = pa14_U_df.T.corr()

# +
# Create a similarity matrix usingn the threshold defined above
# The similarity matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by by the correlation score, otherwise the value is set to 0
# pao1_corr_spell[pao1_corr_spell.abs() < corr_threshold] = 0.0
# pa14_corr_spell[pa14_corr_spell.abs() < corr_threshold] = 0.0

# +
# pao1_corr_spell.head()

# +
# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_corr_spell.abs(), cmap="viridis")
h1.fig.suptitle(f"Correlation of PAO1 genes using threshold={corr_threshold}")

# Save
pao1_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_corr_{corr_threshold}_clustermap.png"
)
# h1.savefig(pao1_clustermap_filename, dpi=300)

# +
plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_corr_spell.abs(), cmap="viridis")
h2.fig.suptitle(f"Correlation of PA14 genes using threshold={corr_threshold}")

# Save
pa14_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_corr_{corr_threshold}_clustermap.png"
)
# h2.savefig(pa14_clustermap_filename, dpi=300)
# -

# ## SEEK

# Correlation of U
pao1_corr_seek = pao1_compendium.corr()
pa14_corr_seek = pa14_compendium.corr()

# +
# Subtract mean correlation score for gene
pao1_corr_mean = pao1_corr_seek.mean()
pa14_corr_mean = pa14_corr_seek.mean()

pao1_corr_mean
# -

pao1_corr_corrected = pao1_corr_seek - pao1_corr_mean
pa14_corr_corrected = pa14_corr_seek - pa14_corr_mean

pao1_corr_seek.head()

# Notice that the corrected correlation matrix is not symmetric anymore because we
# are subtracting the mean for each row
pao1_corr_corrected.head()

pa14_corr_corrected.head()

# +
# Create a similarity matrix usingn the threshold defined above
# The similarity matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by by the correlation score, otherwise the value is set to 0
# pao1_corr_corrected[pao1_corr_corrected.abs() < corr_threshold] = 0.0
# pa14_corr_corrected[pa14_corr_corrected.abs() < corr_threshold] = 0.0

# pao1_corr_corrected.head()

# +
# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_corr_corrected.abs(), cmap="viridis")
h1.fig.suptitle(f"Correlation of PAO1 genes using threshold={corr_threshold}")

# Save
pao1_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_corr_{corr_threshold}_clustermap.png"
)
# h1.savefig(pao1_clustermap_filename, dpi=300)

# +
plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_corr_corrected.abs(), cmap="viridis")
h2.fig.suptitle(f"Correlation of PA14 genes using threshold={corr_threshold}")

# Save
pa14_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_corr_{corr_threshold}_clustermap.png"
)
# h2.savefig(pa14_clustermap_filename, dpi=300)
# -

"""# Save
pao1_corr_filename = f"pao1_corr_{corr_threshold}.tsv"
pa14_corr_filename = f"pa14_corr_{corr_threshold}.tsv"
pao1_corr.to_csv(os.path.join(paths.LOCAL_DATA_DIR, pao1_corr_filename), sep="\t")
pa14_corr.to_csv(os.path.join(paths.LOCAL_DATA_DIR, pa14_corr_filename), sep="\t")"""

# **Takeaway:**
#
# Here we are visualizing the clustering of raw correlation scores where values < `corr_threshold` are set to 0. If we compare the clustermap results in this notebook with [1a_get_network_communities_complex.ipynb](1a_get_network_communities_complex.ipynb) where we cluster the on the Topological Overlap Matrix (TOM) we see:
# * Clustering pattern between using raw correlation score vs TOM is similar. TOM is considering secondary relationships (i.e. gene _i_ and _j_ are similar if they are linked in the adjacency matrix and gene _i_ is connected to all the neighbors of gene _j_)
# * At thresholds 0.5, 0.6 there seems to be 1 large cluster, some very smaller clusters, then all other genes that are below the threshold
# * As we increase the threshold to 0.8 and 0.9, this very large cluster is broken up into more equal sized smaller clusters
# * In terms of distance, for threshold of 0.5, 0.6 high density regions look like those > 20. For thresholds 0.7-0.9, high density regions look like those > 15.
#
# Overall, clustering may make sense using higher thresholds (0.8, 0.9) and excluding the community containing the remaining genes. However then we are not left with many genes, so perhaps it makes sense to consider a different similarity/correlation metric to use?
# * Looking at the pair plots of the raw expression data (estimated counts - an estimate of the number of reads drawn from this transcript given the transcript’s relative abundance and length). There is a tendency for genes to have a long right tail where some genes have a spike at 0 and some do not. In this case Spearman correlation might be more appropriate here.
#
# After meeting with Casey, the main takeaway is that:
# * It appears that the TOM matrix is not as sensitive since its setting nearby genes to have a score of 1.
# Using just the correlation score is grouping genes into one large cluster.
# * To dappen the overwhelming signal of highly correlated genes we will look into applying these two methods:
# https://pubmed.ncbi.nlm.nih.gov/17724061/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4768301/
