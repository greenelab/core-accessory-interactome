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

# # Stability vs operon size
#
# Let's check if co-operonic genes tend to be more stable compared to non co-operonic genes. This will be something we want to note in the manuscript text.
#
# Intuitively, it would make sense that _____

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import textwrap
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import paths, utils, annotations

random.seed(1)
# -

# ## Load data

# Input similarity scores and annotations filenames
# Since the results are similar we only need to look at the scores for one strain type
pao1_similarity_filename = "pao1_core_similarity_associations_final_spell.tsv"
pa14_similarity_filename = "pa14_core_similarity_associations_final_spell.tsv"

# Import df
pao1_similarity = pd.read_csv(pao1_similarity_filename, sep="\t", index_col=0, header=0)
pa14_similarity = pd.read_csv(pa14_similarity_filename, sep="\t", index_col=0, header=0)

print(pao1_similarity.shape)
pao1_similarity.head()

print(pa14_similarity.shape)
pa14_similarity.head()

# ## Load operon infromation

pao1_operon_filename = paths.PAO1_OPERON
pa14_operon_filename = paths.PA14_OPERON

pao1_operon = annotations.load_format_operons(pao1_operon_filename)
pa14_operon = annotations.load_format_operons(pa14_operon_filename)

print(pao1_operon.shape)
pao1_operon.head()

print(pa14_operon.shape)
pa14_operon.head()

# Add operon size column
pao1_operon_size = pao1_operon.value_counts().to_frame("operon_size")
pa14_operon_size = pa14_operon.value_counts().to_frame("operon_size")

# Merge gene id and operon size
pao1_operon_info = pao1_operon.merge(
    pao1_operon_size, left_on="operon_name", right_index=True
)
pa14_operon_info = pa14_operon.merge(
    pa14_operon_size, left_on="operon_name", right_index=True
)
pao1_operon_info.head(10)

pa14_operon_info.head(10)

#  ## Merge operon metadata and stability data

pao1_similarity_operon = pao1_similarity.merge(
    pao1_operon_info, left_index=True, right_index=True, how="left"
)
pa14_similarity_operon = pa14_similarity.merge(
    pa14_operon_info, left_index=True, right_index=True, how="left"
)

print(pao1_similarity_operon.shape)
pao1_similarity_operon.head()

print(pa14_similarity_operon.shape)
pa14_similarity_operon.head()

# ## Plot

# +
# Correlation plot of stability vs operon size
# This will only show the relationship for genes within an operon
# Plot correlation
fig = sns.jointplot(
    data=pao1_similarity_operon,
    x="Transcriptional similarity across strains",
    y="operon_size",
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes")

fig.set_axis_labels(
    "Transcriptional similarity",
    "Operon sizse",
    fontsize=14,
    fontname="Verdana",
)
fig.fig.suptitle(
    "Stability vs operon size (PAO1)", fontsize=16, fontname="Verdana", y=0.9, x=0.45
)

# +
fig = sns.jointplot(
    data=pa14_similarity_operon,
    x="Transcriptional similarity across strains",
    y="operon_size",
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes")

fig.set_axis_labels(
    "Transcriptional similarity",
    "Operon sizse",
    fontsize=14,
    fontname="Verdana",
)
fig.fig.suptitle(
    "Stability vs operon size (PA14)", fontsize=16, fontname="Verdana", y=0.9, x=0.45
)

# +
# Compare the median distribution of transcriptional similarity score for
# genes in operons vs genes not within an operon
pao1_similarity_non_operon = pao1_similarity_operon[
    pao1_similarity_operon["operon_name"].isna()
].index
pao1_similarity_operon_only = pao1_similarity_operon[
    ~pao1_similarity_operon["operon_name"].isna()
].index

pa14_similarity_non_operon = pa14_similarity_operon[
    pa14_similarity_operon["operon_name"].isna()
].index
pa14_similarity_operon_only = pa14_similarity_operon[
    ~pa14_similarity_operon["operon_name"].isna()
].index
# -

pao1_similarity_operon.loc[pao1_similarity_non_operon][
    "Transcriptional similarity across strains"
]

pao1_similarity_operon.loc[pao1_similarity_operon_only][
    "Transcriptional similarity across strains"
]

print(pao1_similarity_non_operon.shape)
print(pao1_similarity_operon_only.shape)

print(pa14_similarity_non_operon.shape)
print(pa14_similarity_operon_only.shape)

# +
# Add label for genes being within an operon or not
pao1_similarity_operon.loc[pao1_similarity_non_operon, "label"] = "non-operon"
pao1_similarity_operon.loc[pao1_similarity_operon_only, "label"] = "operon"

pa14_similarity_operon.loc[pa14_similarity_non_operon, "label"] = "non-operon"
pa14_similarity_operon.loc[pa14_similarity_operon_only, "label"] = "operon"

# +
# Plot coverage distribution given list of generic coverage, specific coverage
pao1_box_fig = sns.boxplot(
    data=pao1_similarity_operon,
    x="label",
    y="Transcriptional similarity across strains",
    notch=True,
    palette=["#81448e", "lightgrey"],
)

pao1_box_fig.set_xlabel(None)
pao1_box_fig.set_xticklabels(
    ["non-operon genes", "operon genes"], fontsize=14, fontname="Verdana"
)
pao1_box_fig.set_ylabel(
    textwrap.fill("Transcriptional similarity", width=30),
    fontsize=14,
    fontname="Verdana",
)
pao1_box_fig.tick_params(labelsize=14)
pao1_box_fig.set_title(
    "Transcriptional similarity if in operon (PAO1)", fontsize=16, fontname="Verdana"
)

# +
pa14_box_fig = sns.boxplot(
    data=pa14_similarity_operon,
    x="label",
    y="Transcriptional similarity across strains",
    notch=True,
    palette=["#81448e", "lightgrey"],
)

pa14_box_fig.set_xlabel(None)
pa14_box_fig.set_xticklabels(
    ["non-operon genes", "operon genes"], fontsize=14, fontname="Verdana"
)
pa14_box_fig.set_ylabel(
    textwrap.fill("Transcriptional similarity", width=30),
    fontsize=14,
    fontname="Verdana",
)
pa14_box_fig.tick_params(labelsize=14)
pa14_box_fig.set_title(
    "Transcriptional similarity if in operon (PA14)", fontsize=16, fontname="Verdana"
)
