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

# # Diversity of strains
#
# To visualize the diversity within strains compared to between strains, we overlay the PCA plot of the samples in each compendium.
#
# Then to quantify the difference we compare the difference in centroid between the PAO1 and PA14 compendium with the spread across the samples within the compendium.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import plotnine as pn
import seaborn as sns
import numpy as np
from scipy.spatial.distance import pdist
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scripts import paths

# +
# Compendia BEFORE processing (normalized filtered counts with all samples mapped)
pao1_all_expression_filename = paths.PAO1_GE
pa14_all_expression_filename = paths.PA14_GE

# Compendia AFTER processing (normalized filtered counts with samples split into strain type)
pao1_strain_expression_filename = paths.PAO1_COMPENDIUM
pa14_strain_expression_filename = paths.PA14_COMPENDIUM

# File containing table to map sample id to strain name
# sample_to_strain_filename = paths.SAMPLE_TO_STRAIN

# +
# Load expression data
pao1_all_expression = pd.read_csv(
    pao1_all_expression_filename, index_col=0, header=0, sep="\t"
)
pa14_all_expression = pd.read_csv(
    pa14_all_expression_filename, index_col=0, header=0, sep="\t"
)

pao1_strain_expression = pd.read_csv(
    pao1_strain_expression_filename, index_col=0, header=0, sep="\t"
)
pa14_strain_expression = pd.read_csv(
    pa14_strain_expression_filename, index_col=0, header=0, sep="\t"
)
# -

# Confirm that BEFORE compendia have 2,333 samples
print(pao1_all_expression.shape)
pao1_all_expression.head()

print(pa14_all_expression.shape)
pa14_all_expression.head()

# ## Get labels for which samples are in PAO1 and PA14 compendia

# Confirm that PAO1 has 890 samples and PA14 has 505 samples
print(pao1_strain_expression.shape)
print(pa14_strain_expression.shape)

# +
# Add labels for samples that belong to compendium
pao1_strain_expression["compendium"] = "PAO1"
pao1_strain_label = pao1_strain_expression["compendium"].to_frame()

print(pao1_strain_label.shape)
pao1_strain_label.head()

# +
pa14_strain_expression["compendium"] = "PA14"
pa14_strain_label = pa14_strain_expression["compendium"].to_frame()

print(pa14_strain_label.shape)
pa14_strain_label.head()

# +
# Merge label dfs
strain_labels = pd.concat([pao1_strain_label, pa14_strain_label])

assert strain_labels.shape[0] == pao1_strain_label.shape[0] + pa14_strain_label.shape[0]
print(strain_labels.shape)
strain_labels.head()
# -

# ### Format expression compendia
#
# Train PCA on PAO1 compendium and apply to PA14 compendium

# ### Low dimensional embedding

# +
# Try 0-1 scaling for PCA
# scaler = MinMaxScaler()
scaler = StandardScaler()

# Fitting
normalized_pao1_expression = scaler.fit_transform(pao1_all_expression)
normalized_pa14_expression = scaler.fit_transform(pa14_all_expression)

normalized_pao1_expression_df = pd.DataFrame(
    normalized_pao1_expression,
    columns=pao1_all_expression.columns,
    index=pao1_all_expression.index,
)

normalized_pa14_expression_df = pd.DataFrame(
    normalized_pa14_expression,
    columns=pa14_all_expression.columns,
    index=pa14_all_expression.index,
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
    index=pao1_all_expression.index,
    columns=["1", "2"],
)

# Format
pa14_pca_encoded_df = pd.DataFrame(
    data=pa14_pca_encoded,
    index=pa14_all_expression.index,
    columns=["1", "2"],
)

# Add strain labels
pao1_pca_encoded_label = pao1_pca_encoded_df.merge(
    strain_labels, left_index=True, right_index=True
)

print(pao1_pca_encoded_label.shape)
pao1_pca_encoded_label.head()

pa14_pca_encoded_label = pa14_pca_encoded_df.merge(
    strain_labels, left_index=True, right_index=True
)

# ### Plot

# Colors
edge_colors = {
    # "Clinical Isolate": "#89A45E",
    "PA14": "#895881",
    # "PAK": "#EF8B46",
    "PAO1": "#C6A9B5",
    # "NA": "#D8DAEB",
    # "Other": "#808080"
}

# +
# Plot gene expression in PAO1 reference
fig1 = pn.ggplot(pao1_pca_encoded_label, pn.aes(x="1", y="2"))
fig1 += pn.geom_point(pn.aes(color="compendium"), alpha=0.3, size=3, stroke=0.8)
fig1 += pn.scale_color_manual(values=edge_colors)
fig1 += pn.labs(
    x="PCA 1",
    y="PCA 2",
    title="Expression using PAO1 reference",
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

fig1.save("pa_pao1_ref_pca.svg", format="svg", dpi=300)

# +
# Plot gene expression in PA14 reference
fig3 = pn.ggplot(pa14_pca_encoded_label, pn.aes(x="1", y="2"))
fig3 += pn.geom_point(pn.aes(color="compendium"), alpha=0.3, size=3, stroke=0.8)
fig3 += pn.scale_color_manual(values=edge_colors)
fig3 += pn.labs(
    x="PCA 1",
    y="PCA 2",
    title="Expression using PA14 reference",
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

fig3.save("pa_pa14_ref_pca.svg", format="svg", dpi=300)
# -

# ## Calculate similarity between compendium

# +
# What is the variance explained
pca_pao1 = PCA(n_components=200)
encoded_pca_pao1 = pca_pao1.fit_transform(normalized_pao1_expression_df)

exp_var_pca_pao1 = pca_pao1.explained_variance_ratio_
cum_sum_eigenvalues_pao1 = np.cumsum(exp_var_pca_pao1)

plt.step(
    range(0, len(cum_sum_eigenvalues_pao1)),
    cum_sum_eigenvalues_pao1,
    where="mid",
    label="Cumulative explained variance",
)
plt.ylabel("Explained variance ratio")

# +
pca_pa14 = PCA(n_components=200)
encoded_pca_pa14 = pca_pa14.fit_transform(normalized_pa14_expression_df)

exp_var_pca_pa14 = pca_pa14.explained_variance_ratio_
cum_sum_eigenvalues_pa14 = np.cumsum(exp_var_pca_pa14)

plt.step(
    range(0, len(cum_sum_eigenvalues_pa14)),
    cum_sum_eigenvalues_pa14,
    where="mid",
    label="Cumulative explained variance",
)
plt.ylabel("Explained variance ratio")
# -

# Looks like 200 PCs explains 90% of the variance. So let's use these top 200 PCs to calculate the difference in centroid

# Get centroid
pao1_centroid = encoded_pca_pao1.mean(axis=0)
pa14_centroid = encoded_pca_pa14.mean(axis=0)

pao1_centroid = pao1_centroid.reshape(-1, 1)
pa14_centroid = pa14_centroid.reshape(-1, 1)

both_centroid = np.vstack((pao1_centroid.T, pa14_centroid.T))

# +
# Calculate distance between centroids
mean_dist = pdist(both_centroid)

mean_dist
# -

# Variance of PAO1 compendium
pca_pao1.explained_variance_.sum()

# Variance of PA14 compendium
pca_pa14.explained_variance_.sum()

encoded_pca_pa14.var(axis=0)
