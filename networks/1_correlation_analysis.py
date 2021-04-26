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
# This notebook performs correlation analysis to compare the similarity between genes and applies different threshold cutoffs to determine the strength of connection between genes

# %load_ext autoreload
# %autoreload 2
import os
import pandas as pd
import plotnine as pn
import seaborn as sns
import matplotlib.pyplot as plt
import umap
from scipy.spatial.distance import pdist, squareform
from core_acc_modules import paths

# ## Set user parameters
#
# For now we will vary the correlation threshold (`corr_threshold`) but keep the other parameters consistent
#
# We will run this notebook for each threshold parameter

# +
# Params
corr_threshold = 0.9

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

# ## Get similarity between genes
#
# To determine if genes are similar, we will calculate the correlation between genes and apply our threshold cutoff. When we apply our threshold, any scores that are below the threshold will be set to 0 (i.e. using a threshold of 0.5 means that gene pairs that have correlation scores of 0.52 and 0.78 will be left as is but a gene pair with a correlation score of 0.48 will be set to 0).

# Get perason correlation
# This correlation matrix will represent the concordance
# between two gene expression profiles
pao1_corr = pao1_compendium.corr()
pa14_corr = pa14_compendium.corr()

pao1_corr.head()

# +
# Create a similarity matrix usingn the threshold defined above
# The similarity matrix will determine the strength of the connection between two genes
# If the concordance is strong enough (i.e. above the threshold), then
# the genes are connected by by the correlation score, otherwise the value is set to 0
pao1_corr[pao1_corr.abs() < corr_threshold] = 0.0
pa14_corr[pa14_corr.abs() < corr_threshold] = 0.0

pao1_corr.head()
# -

pa14_corr.head()

# ## Plot distribution of pairwise distances
#
# This will particularly help to inform the parameters we use for DBSCAN, which is density based

# Get distribution of pairwise distances to determine a cutoff defining what a dense region should be
f1 = sns.displot(pdist(pao1_corr))
plt.title("Distribution of pairwise distances for PAO1 genes")

f2 = sns.displot(pdist(pa14_corr))
plt.title("Distribution of pairwise distances for PA14 genes")

# ## Plot correlation
#
# We will plot a heatmap and umap of the correlations to identify clusters, which should help to inform the parameters for hierarchal clustering - i.e. how many clusters can we expect?

# +
# Plot heatmap
plt.figure(figsize=(20, 20))
h1 = sns.clustermap(pao1_corr.abs(), cmap="viridis")
h1.fig.suptitle(f"Correlation of PAO1 genes using threshold={corr_threshold}")

# Save
pao1_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_corr_{corr_threshold}_clustermap.png"
)
h1.savefig(pao1_clustermap_filename, dpi=300)

# +
plt.figure(figsize=(20, 20))
h2 = sns.clustermap(pa14_corr.abs(), cmap="viridis")
h2.fig.suptitle(f"Correlation of PA14 genes using threshold={corr_threshold}")

# Save
pa14_clustermap_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_corr_{corr_threshold}_clustermap.png"
)
h2.savefig(pa14_clustermap_filename, dpi=300)

# +
model_pao1 = umap.UMAP(random_state=123).fit(pao1_corr)

pao1_encoded = model_pao1.transform(pao1_corr)

pao1_encoded_df = pd.DataFrame(
    data=pao1_encoded,
    index=pao1_corr.index,
    columns=["1", "2"],
)

# +
model_pa14 = umap.UMAP(random_state=123).fit(pa14_corr)

pa14_encoded = model_pa14.transform(pa14_corr)

pa14_encoded_df = pd.DataFrame(
    data=pa14_encoded,
    index=pa14_corr.index,
    columns=["1", "2"],
)

# +
# Plot PAO1
u1 = pn.ggplot(pao1_encoded_df, pn.aes(x="1", y="2"))
u1 += pn.geom_point(pn.aes(alpha=0.1))
u1 += pn.labs(x="UMAP 1", y="UMAP 2", title="Correlation of PAO1 genes")
u1 += pn.theme_bw()
u1 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)

print(u1)

# +
# Plot PA14
u2 = pn.ggplot(pa14_encoded_df, pn.aes(x="1", y="2"))
u2 += pn.geom_point(pn.aes(alpha=0.2))
u2 += pn.labs(x="UMAP 1", y="UMAP 2", title="Correlation of PA14 genes")
u2 += pn.theme_bw()
u2 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)

print(u2)
# -

# Save
pao1_corr_filename = f"pao1_corr_{corr_threshold}.tsv"
pa14_corr_filename = f"pa14_corr_{corr_threshold}.tsv"
pao1_corr.to_csv(os.path.join(paths.LOCAL_DATA_DIR, pao1_corr_filename), sep="\t")
pa14_corr.to_csv(os.path.join(paths.LOCAL_DATA_DIR, pa14_corr_filename), sep="\t")

# **Observations:**
# * At a threshold of 0.5, there is 1 very large cluster and a few very small clusters and then everything else for PAO1 genes. For PA14 genes, there is 1 large cluster and everything else. High density regions look like those > 20.
# * At a threshold of 0.6, there is 1 very large cluster and everything else for PAO1 genes. For PA14 genes there looks to be 1 large cluster and 1 smaller cluster and then everything else. High density regions look like those > 20.
# * At a threshold of 0.7, High density regions look like those > 15.
# * At a threshold of 0.8, High density regions look like those > 15.
# * At a threshold of 0.9, High density regions look like those > 15.
#
# Overall there are smaller clusters at a higher threshold, which we would expect since we wouldn't expect as many genes to have such a high correlation. We will use the number of clusters observed to help validate the clustering results in the next notebook.
#
# What do we do about the remaining group of genes (i.e. the ones with 0 correlation)?
