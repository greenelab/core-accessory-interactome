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

# # Transformation + correlation analysis
#
# This notebook examines the correlation structure in the gene expression data generated in [1_create_compendia.ipynb](../processing/1_create_compendia.ipynb).
#
# When we performed clustering on the correlation matrices (using Pearson correlation) we found that pairs of genes had either very high correlation scores (>0.5) or very low correlation scores (<0.1). As a result gene pairs that were highly correlated clustered into a single large module. This clustering pattern is not ideal for a couple of reasons:
# 1. Given that our goal is to examine the composition of gene groups, having all genes within one module does not allow us to do this
# 2. These highly correlated group of genes are likely masking other relevant specific signals/relationships in the data
#
# Here we will apply transformations of the raw gene expression data before applying clustering to see if this helps to correct for this high level of correlation

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
from sklearn import preprocessing
from core_acc_modules import paths

# ## Set user parameters
#
# For now we will vary the correlation threshold (`corr_threshold`) but keep the other parameters consistent
#
# We will run this notebook for each threshold parameter

# +
# Params -- REMOVE???
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

# ## Scaling + correlation
#
# Try different processing of the raw gene expression data and then applying correlation:
# * log transform + correlation
# * Normalize + correlation

# ### Log transform + correlation

# log transform data
pao1_compendium_log10 = np.log10(pao1_compendium)
pa14_compendium_log10 = np.log10(pa14_compendium)

# Correlation
pao1_corr_log10 = pao1_compendium_log10.corr()
pa14_corr_log10 = pa14_compendium_log10.corr()

# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_corr_log10.abs(), cmap="viridis")
h1.fig.suptitle("Correlation of log10 transformed PAO1 genes", y=1.05)

plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_corr_log10.abs(), cmap="viridis")
h2.fig.suptitle("Correlation of log10 transformed PA14 genes", y=1.05)

# ### 0-1 normalize + correlation

# +
# 0-1 normalize per gene
scaler = preprocessing.MinMaxScaler()

pao1_scaled = scaler.fit_transform(pao1_compendium)
pao1_scaled_df = pd.DataFrame(
    pao1_scaled, columns=pao1_compendium.columns, index=pao1_compendium.index
)

scaler2 = preprocessing.MinMaxScaler()
pa14_scaled = scaler2.fit_transform(pa14_compendium)
pa14_scaled_df = pd.DataFrame(
    pa14_scaled, columns=pa14_compendium.columns, index=pa14_compendium.index
)
# -

sns.displot(pao1_scaled_df["PA0001"])

# Correlation
pao1_corr_normalized = pao1_scaled_df.corr()
pa14_corr_normalized = pa14_scaled_df.corr()

# Plot heatmap
plt.figure(figsize=(20, 20))
h3 = sns.clustermap(pao1_corr_normalized.abs(), cmap="viridis")
h3.fig.suptitle("Correlation of 0-1 normalized PAO1 genes", y=1.05)

plt.figure(figsize=(20, 20))
h4 = sns.clustermap(pa14_corr_normalized.abs(), cmap="viridis")
h4.fig.suptitle("Correlation of 0-1 normalized PA14 genes", y=1.05)

# **Takeaway:**
#
# * Using log transform, it looks like there are some smaller modules forming but still one very large module for PAO1. For the pA14 data, it looks like there are more equal-sized modules
# * Normalizing the data, it looks like there is predominantly one large module for both PAO1 and PA14.
