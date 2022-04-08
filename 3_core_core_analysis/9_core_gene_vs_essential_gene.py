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

# # Core gene vs essential gene
#
# [Poulsen et al.](https://www.pnas.org/doi/10.1073/pnas.1900570116) published a study identifying and characterizing essential (i.e. critical for the organism's survival) core genes. Here we want to compare the 321 essential genes identified compared to our most stable core set of genes.

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
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from scripts import paths, utils, annotations

random.seed(1)
# -

# ## Load data

# +
# Input similarity scores
our_similarity_filename = "all_core_similarity_associations_final_spell.tsv"

# Downloaded supplemental table S6
published_essential_filename = "pnas_1900570116_sd06.csv"
# -

# Import df
our_similarity = pd.read_csv(our_similarity_filename, sep="\t", index_col=0, header=0)
published_essential = pd.read_csv(published_essential_filename, index_col=0, header=0)

print(our_similarity.shape)
our_similarity.head()

print(published_essential.shape)
published_essential.head()

# +
# Essential genes are only those with an E in the "This study"
published_essential = published_essential[published_essential["This study"] == "E"]

# Check that there are 321 genes
print(published_essential.shape)
published_essential.head()
# -

# ## Compare most stable core genes vs core essential genes

# Get most stable core genes
most_stable_core = our_similarity.loc[
    our_similarity["label"] == "most stable", "PA14 homolog id"
]

print(len(most_stable_core))
most_stable_core

# Get core essential genes
pub_core_essential = published_essential["PA14_ID (for reference)"]

# +
# Compare
most_stable_venn = venn2(
    [set(most_stable_core), set(pub_core_essential)],
    set_labels=("Most stable core", "Essential core"),
)

most_stable_venn.get_patch_by_id("11").set_color("purple")
most_stable_venn.get_patch_by_id("11").set_edgecolor("none")
most_stable_venn.get_patch_by_id("11").set_alpha(0.3)
most_stable_venn.get_patch_by_id("01").set_color("blue")
most_stable_venn.get_patch_by_id("01").set_edgecolor("none")
most_stable_venn.get_patch_by_id("01").set_alpha(0.3)

plt.title("Most stable core vs core essental", fontsize=16, fontname="Verdana")
for text in most_stable_venn.set_labels:
    text.set_fontsize(14)
    text.set_fontname("Verdana")

for text in most_stable_venn.subset_labels:
    text.set_fontsize(12)
    text.set_fontname("Verdana")
# -

# % of core essential genes that are most stable
109 / 267

# ### Explore genes a little

our_similarity = our_similarity.set_index("PA14 homolog id")

# +
# Genes that are shared
shared_genes = set(most_stable_core).intersection(set(pub_core_essential))

our_similarity.loc[shared_genes]

# +
# Genes that are only most stable
only_stable_genes = set(most_stable_core).difference(set(pub_core_essential))

our_similarity.loc[only_stable_genes]

# +
# Genes that are only essential
only_essential_genes = set(pub_core_essential).difference(set(most_stable_core))
shared_only_essential_genes = only_essential_genes.intersection(our_similarity.index)

our_similarity.loc[shared_only_essential_genes]
