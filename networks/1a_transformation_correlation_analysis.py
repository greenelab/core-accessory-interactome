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

# # Scaling & correlation
#
# Try different processing of the raw gene expression data and then applying correlation:
# * log transform then correlation
# * Normalize then correlation
#
# ## update text based on papers
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

# Plot heatmap
plt.figure(figsize=(20, 20))
o1 = sns.clustermap(pao1_corr_original.abs(), cmap="viridis")
o1.fig.suptitle("Correlation of raw PAO1 genes", y=1.05)

# Plot heatmap
plt.figure(figsize=(20, 20))
o2 = sns.clustermap(pa14_corr_original.abs(), cmap="viridis")
o2.fig.suptitle("Correlation of raw PA14 genes", y=1.05)

# ## Log transform + correlation

# log transform data
# Note: add 1 to avoid -inf and so 0 corresponds to those with 0 counts
pao1_compendium_log10 = np.log10(1 + pao1_compendium)
pa14_compendium_log10 = np.log10(1 + pa14_compendium)

sns.displot(pao1_compendium_log10["PA0001"])

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

# ## 0-1 normalize + correlation

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
# Here we are using the Pearson correlation metric, which is defined by the formula:
#
# $$ \frac{\sigma_{XY}}{\sigma_X \sigma_Y}$$
#
# In other words, the correlation is calculated based on the covariance between X,Y (i.e. how X jointly varies with Y) and the scaled by the standard deviation of X and Y (i.e. how variable X and Y are). So this Pearson correlation is testing for linear associations between X and Y.
#
# * We find the same correlation patterns if we compare the raw data vs normalized data. This is because the normalization we perform is scaling the values linearly.
# * By comparison, if we apply a non-linear scaling like log transform, the correlation values are changed. In fact, the correlation patterns using log data reveals more smaller modules compared to using raw data. After applying the log transform to our data, the data will be compressed (i.e. values in the right tail will be scaled down and these high values are what are creating the very dominant high correlation scores. Looking at the equation above the correlation score will be high if the covariance is high. And the covariance will be high if gene X and Y have long right tails and low means which they tend to as seen in the displots above).
