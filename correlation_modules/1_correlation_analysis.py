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
from core_acc_modules import paths, utils

# ## Set user parameters
#
# Here we set the number of PCs or singular vectors to use. We are starting with 300 since this is what [eADAGE](https://pubmed.ncbi.nlm.nih.gov/28711280/) used.

# +
# Params
num_singular_values = 300
num_singular_values_log = 100

# Which subset of genes to consider: core, acc, all
subset_genes = "all"
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

# ## Correlation of raw gene expression data
#
# Here is the correlation of the raw data without any malnipulations. This will serve as a reference to compare the correlations below where applied corrections to the correlations to account for the dominant signal described above.

# Correlation
pao1_corr_original = pao1_compendium.corr()
pa14_corr_original = pa14_compendium.corr()

pao1_corr_original.head()

pa14_corr_original.head()

# +
# %%time
# Plot heatmap
o1 = sns.clustermap(pao1_corr_original.abs(), cmap="viridis", figsize=(20, 20))
o1.fig.suptitle("Correlation of raw PAO1 genes", y=1.05)

# Save
pao1_pearson_filename = os.path.join(
    paths.LOCAL_DATA_DIR, "pao1_pearson_clustermap.png"
)
o1.savefig(pao1_pearson_filename, dpi=300)

# +
# Plot heatmap
o2 = sns.clustermap(pa14_corr_original.abs(), cmap="viridis", figsize=(20, 20))
o2.fig.suptitle("Correlation of raw PA14 genes", y=1.05)

# Save
pa14_pearson_filename = os.path.join(
    paths.LOCAL_DATA_DIR, "pa14_pearson_clustermap.png"
)
o2.savefig(pa14_pearson_filename, dpi=300)
# -

# Save original correlation matrices
pao1_pearson_mat_filename = os.path.join(paths.LOCAL_DATA_DIR, "pao1_pearson_mat.tsv")
pa14_pearson_mat_filename = os.path.join(paths.LOCAL_DATA_DIR, "pa14_pearson_mat.tsv")
pao1_corr_original.to_csv(pao1_pearson_mat_filename, sep="\t")
pa14_corr_original.to_csv(pa14_pearson_mat_filename, sep="\t")

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

print(pao1_compendium_T.shape)
print(pao1_U.shape, pao1_s.shape, pao1_Vh.shape)

print(pa14_compendium_T.shape)
print(pa14_U.shape, pa14_s.shape, pa14_Vh.shape)

# Convert ndarray to df to use corr()
pao1_U_df = pd.DataFrame(data=pao1_U, index=pao1_compendium_T.index)
pa14_U_df = pd.DataFrame(data=pa14_U, index=pa14_compendium_T.index)

# +
# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U

pao1_corr_log_spell = pao1_U_df.iloc[:, :num_singular_values_log].T.corr()
pa14_corr_log_spell = pa14_U_df.iloc[:, :num_singular_values_log].T.corr()
# -

# Select subset of genes
if subset_genes == "core":
    pao1_corr_log_spell = pao1_corr_log_spell.loc[pao1_core, pao1_core]
    pa14_corr_log_spell = pa14_corr_log_spell.loc[pa14_core, pa14_core]
elif subset_genes == "acc":
    pao1_corr_log_spell = pao1_corr_log_spell.loc[pao1_acc, pao1_acc]
    pa14_corr_log_spell = pa14_corr_log_spell.loc[pa14_acc, pa14_acc]

# +
# Plot heatmap
h1a = sns.clustermap(pao1_corr_log_spell.abs(), cmap="viridis", figsize=(20, 20))
h1a.fig.suptitle(
    f"log transform + SPELL corrected using {num_singular_values_log} vectors (PAO1)",
    y=1.05,
)

# Save
pao1_log_spell_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_log_spell_clustermap.png"
)
h1a.savefig(pao1_log_spell_filename, dpi=300)

# +
h2a = sns.clustermap(pa14_corr_log_spell.abs(), cmap="viridis", figsize=(20, 20))
h2a.fig.suptitle(
    f"log transformed + SPELL corrected using {num_singular_values_log} vectors (PA14)",
    y=1.05,
)

# Save
pa14_log_spell_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_{subset_genes}_log_spell_clustermap.png"
)
h2a.savefig(pa14_log_spell_filename, dpi=300)
# -

# ## Plot distribution of pairwise distances
#
# This will particularly help to inform the parameters we use for DBSCAN, which is density based. Here we looking at the distribution of both global distances and local distances. Global distances are defined using `pdist`, which takes the pairwise Euclidean distance of each of the correlation vectors (so the distance between gene `p` and gene `q` is based on the difference in correlation between `p` and all other genes, and `q` and all other genes). Whereas the local distance is defined as 1 - |correlation(`p`, `q`)|

# Get distribution of pairwise distances to determine a cutoff defining what a dense region should be
f1 = sns.displot(pdist(pao1_corr_log_spell))
plt.title("Distribution of pairwise distances for PAO1 genes")

f2 = sns.displot(pdist(pa14_corr_log_spell))
plt.title("Distribution of pairwise distances for PA14 genes")

# +
pao1_local_dist = 1 - pao1_corr_log_spell.abs()
pao1_local_dist = pao1_local_dist.where(
    np.triu(np.ones(pao1_local_dist.shape), k=1).astype(np.bool)
)
pao1_local_dist = pao1_local_dist.stack().reset_index()
pao1_local_dist.columns = ["Row", "Column", "Value"]

pao1_local_dist.head(10)
# -

f3 = sns.displot(pao1_local_dist["Value"])
plt.title("Distribution of pairwise distances for PAO1 genes")

# +
pa14_local_dist = 1 - pa14_corr_log_spell.abs()
pa14_local_dist = pa14_local_dist.where(
    np.triu(np.ones(pa14_local_dist.shape), k=1).astype(np.bool)
)
pa14_local_dist = pa14_local_dist.stack().reset_index()
pa14_local_dist.columns = ["Row", "Column", "Value"]

pa14_local_dist.head(10)
# -

f4 = sns.displot(pa14_local_dist["Value"])
plt.title("Distribution of pairwise distances for PA14 genes")

# Save log transform + SPELL correlation matrices
pao1_log_spell_mat_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_{subset_genes}_log_spell_mat.tsv"
)
pa14_log_spell_mat_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_{subset_genes}_log_spell_mat.tsv"
)
pao1_corr_log_spell.to_csv(pao1_log_spell_mat_filename, sep="\t")
pa14_corr_log_spell.to_csv(pa14_log_spell_mat_filename, sep="\t")