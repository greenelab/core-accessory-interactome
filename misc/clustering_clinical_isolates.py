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

# # Clustering of clinical isolates
#
# This notebook examines how clinical isolate samples cluster with PAO1 and PA14 strains using least stable core genes only

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import plotnine as pn
import seaborn as sns
from textwrap import fill
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scripts import paths

# ## Load data

# +
# Expression data pre-binning
pao1_expression_prebin_filename = paths.PAO1_PREBIN_COMPENDIUM
pa14_expression_prebin_filename = paths.PA14_PREBIN_COMPENDIUM

# SRA annotations
sra_annotation_filename = paths.SAMPLE_TO_STRAIN_PROCESSED

# +
# Load SRA labels
sample_to_strain_df = pd.read_csv(
    sra_annotation_filename, sep="\t", index_col=0, header=0
)

sample_to_strain_df.head()
# -

# Load compendium with all 2,333 samples mapped to PAO1, PA14 references
pao1_compendium = pd.read_csv(
    pao1_expression_prebin_filename, sep="\t", index_col=0, header=0
)
pa14_compendium = pd.read_csv(
    pa14_expression_prebin_filename, sep="\t", index_col=0, header=0
)

print(pao1_compendium.shape)
pao1_compendium.head()

# ## Select least stable core genes

# +
# Load transcriptional similarity df
# These are the subset of genes that we will consider
pao1_similarity_scores_filename = (
    "../3_core_core_analysis/pao1_core_similarity_associations_final_spell.tsv"
)
pa14_similarity_scores_filename = (
    "../3_core_core_analysis/pa14_core_similarity_associations_final_spell.tsv"
)

pao1_similarity_scores = pd.read_csv(
    pao1_similarity_scores_filename, sep="\t", header=0, index_col=0
)
pa14_similarity_scores = pd.read_csv(
    pa14_similarity_scores_filename, sep="\t", header=0, index_col=0
)
# -

# Get least stable core genes
pao1_least_stable_genes = list(
    pao1_similarity_scores[pao1_similarity_scores["label"] == "least stable"].index
)
pa14_least_stable_genes = list(
    pa14_similarity_scores[pa14_similarity_scores["label"] == "least stable"].index
)

# Select least stable core genes
pao1_compendium_least = pao1_compendium[pao1_least_stable_genes]
pa14_compendium_least = pa14_compendium[pa14_least_stable_genes]

print(pao1_compendium_least.shape)
pao1_compendium_least.head()

# ## Clustering

# +
# Try 0-1 scaling for PCA
scaler = StandardScaler()

# Fitting
normalized_pao1_expression = scaler.fit_transform(pao1_compendium_least)
normalized_pa14_expression = scaler.fit_transform(pa14_compendium_least)

normalized_pao1_expression_df = pd.DataFrame(
    normalized_pao1_expression,
    columns=pao1_compendium_least.columns,
    index=pao1_compendium_least.index,
)

normalized_pa14_expression_df = pd.DataFrame(
    normalized_pa14_expression,
    columns=pa14_compendium_least.columns,
    index=pa14_compendium_least.index,
)

# +
pca = PCA(n_components=2)
model_pca_pao1 = pca.fit(normalized_pao1_expression_df)

pao1_pca_encoded = model_pca_pao1.transform(normalized_pao1_expression_df)

# +
pca = PCA(n_components=2)
model_pca_pa14 = pca.fit(normalized_pa14_expression_df)

pa14_pca_encoded = model_pca_pa14.transform(normalized_pa14_expression_df)
# -

# Format
pao1_pca_encoded_df = pd.DataFrame(
    data=pao1_pca_encoded,
    index=pao1_compendium_least.index,
    columns=["1", "2"],
)

# Format
pa14_pca_encoded_df = pd.DataFrame(
    data=pa14_pca_encoded,
    index=pa14_compendium_least.index,
    columns=["1", "2"],
)

# Add strain labels
pao1_pca_encoded_label = pao1_pca_encoded_df.merge(
    sample_to_strain_df, left_index=True, right_index=True
)

pa14_pca_encoded_label = pa14_pca_encoded_df.merge(
    sample_to_strain_df, left_index=True, right_index=True
)

print(pao1_pca_encoded_label.shape)
pao1_pca_encoded_label.head()

pao1_pca_encoded_label["Strain type"].value_counts()

# Colors
edge_colors = {
    "Clinical Isolate": "#89A45E",
    "PA14": "#895881",
    # "PAK": "#EF8B46",
    "PAO1": "#C6A9B5",
    # "NA": "#D8DAEB",
    # "Other": "#808080"
}

# +
# Plot gene expression in PAO1 reference

# Only plot PAO1, PA14 and clinical samples
pao1_pca_encoded_label_subset = pao1_pca_encoded_label[
    (pao1_pca_encoded_label["Strain type"] == "PAO1")
    | (pao1_pca_encoded_label["Strain type"] == "PA14")
    | (pao1_pca_encoded_label["Strain type"] == "Clinical Isolate")
]

fig1 = pn.ggplot(pao1_pca_encoded_label_subset, pn.aes(x="1", y="2"))
fig1 += pn.geom_point(pn.aes(color="Strain type"), alpha=0.2, size=1, stroke=0.8)
fig1 += pn.scale_color_manual(values=edge_colors)
fig1 += pn.labs(
    x="PCA 1",
    y="PCA 2",
    title="Least stable core gene expression using PAO1 reference",
)
fig1 += pn.theme_bw()
fig1 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)
fig1 += pn.guides(colour=pn.guide_legend(override_aes={"alpha": 1}))

print(fig1)

# +
# Plot gene expression in PA14 reference

# Only plot PAO1, PA14 and clinical samples
pa14_pca_encoded_label_subset = pao1_pca_encoded_label[
    (pa14_pca_encoded_label["Strain type"] == "PAO1")
    | (pa14_pca_encoded_label["Strain type"] == "PA14")
    | (pa14_pca_encoded_label["Strain type"] == "Clinical Isolate")
]

fig3 = pn.ggplot(pa14_pca_encoded_label_subset, pn.aes(x="1", y="2"))
fig3 += pn.geom_point(pn.aes(color="Strain type"), alpha=0.3, size=1, stroke=0.8)
fig3 += pn.scale_color_manual(values=edge_colors)
fig3 += pn.labs(
    x="PCA 1",
    y="PCA 2",
    title="Least stable core gene expression using PA14 reference",
)
fig3 += pn.theme_bw()
fig3 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)
fig3 += pn.guides(colour=pn.guide_legend(override_aes={"alpha": 1}))
fig3 += pn.guides(fill=pn.guide_legend(override_aes={"alpha": 1}))

print(fig3)
