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

# # Explore module relationships
#
# There are some accessory-accessory modules that we'd like to explore how they are related. To do this, we will plot a heatmap to show how the genes in accessory module A and accessory module B are related.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from core_acc_modules import utils, paths

np.random.seed(1)

# +
# Clustering method
method_name = "affinity"

# Gene subset
gene_subset = "acc"

# Select 2 modules
select_modules = [10, 14]
# -

# ### Load correlation matrix

# Load correlation matrix
pao1_corr_filename = paths.PAO1_CORR_LOG_SPELL_ACC

pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)

# ### Load module membership

pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method_name}_{gene_subset}.tsv"
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", header=0, index_col=0)

select_modules_df = pao1_membership[
    (pao1_membership["module id"] == 10) | (pao1_membership["module id"] == 14)
]

select_modules_df.head()

# ### Heatmap

select_gene_ids = list(select_modules_df.index)

select_corr = pao1_corr.loc[select_gene_ids, select_gene_ids]

select_corr.head()

# Make color labels for module id
my_palette = dict(zip(select_modules_df["module id"].unique(), ["red", "blue"]))
row_colors = select_modules_df["module id"].map(my_palette)

my_palette

# %%time
f = sns.clustermap(
    select_corr.abs(), row_colors=row_colors, cmap="viridis", figsize=(20, 20)
)
f.ax_heatmap.set_xticklabels(f.ax_heatmap.get_xmajorticklabels(), fontsize=16)
f.fig.suptitle(f"Correlation of modules {select_modules}", y=1.05)

# These two accessory modules were of interest to our collaborators. She was interested in how they were related. To get at this, I plotted the heatmap of the correlation matrix, subsetting only for those gene ids that belong to these two modules of interest.
#
# Based on this heatmap, the genes within each module are more closely co-expressed with each other, as expected. There doesn't look to be too many genes that cross over.
