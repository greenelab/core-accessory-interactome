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

# Get distribution of pairwise distances to determine what a dense region should be
squareform(pdist(pao1_corr))

squareform(pdist(pa14_corr))

# ## Plot correlation
#
# We will plot a heatmap and umap of the correlations to identify clusters, which should help to inform the parameters for hierarchal clustering

# +
# Plot heatmap

# +
model_pao1 = umap.UMAP(random_state=123).fit(pao1_corr)

pao1_encoded = model_pao1.transform(pao1_corr)

pao1_encoded_df = pd.DataFrame(
    data=pao1_encoded,
    index=pao1_corr.index,
    columns=["1", "2"],
)

# +
# model_pao1 = pca.fit(normalized_pao1_expression_numeric_df)
model_pa14 = umap.UMAP(random_state=123).fit(pao1_corr)

pao1_encoded = model_pao1.transform(pao1_corr)

pao1_encoded_df = pd.DataFrame(
    data=pao1_encoded,
    index=pao1_corr.index,
    columns=["1", "2"],
)
# -

pao1_encoded_df.head()

# +
# Plot PAO1
fig = pn.ggplot(pao1_encoded_df, pn.aes(x="1", y="2"))
fig += pn.geom_point(pn.aes(alpha=0.3))
fig += pn.labs(x="UMAP 1", y="UMAP 2", title="Correlation of PAO1 genes")
fig += pn.theme_bw()
fig += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)

print(fig)
# -

# Save
pao1_corr_filename = f"pao1_corr_{corr_threshold}.tsv"
pa14_corr_filename = f"pa14_corr_{corr_threshold}.tsv"
pao1_corr.to_csv(os.path.join(paths.LOCAL_DATA_DIR, pao1_corr_filename), sep="\t")
pa14_corr.to_csv(os.path.join(paths.LOCAL_DATA_DIR, pa14_corr_filename), sep="\t")
