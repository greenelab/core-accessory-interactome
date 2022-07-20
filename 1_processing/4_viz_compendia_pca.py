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
# -

# Get variance explain
model_pca_pao1.explained_variance_ratio_

# +
pca = PCA(n_components=2)
model_pca_pa14 = pca.fit(normalized_pa14_expression_df)

pa14_pca_encoded = model_pca_pa14.transform(normalized_pa14_expression_df)
# -

# Get variance explain
model_pca_pa14.explained_variance_ratio_

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

# Plot pao1 and pa14 compendia separately to better visualize the overlap
fig1 = pn.ggplot(
    pao1_pca_encoded_label[pao1_pca_encoded_label["compendium"] == "PAO1"],
    pn.aes(x="1", y="2"),
)
fig1 += pn.geom_point(color="#C6A9B5", alpha=0.3, size=4, stroke=0.8)
fig1 += pn.labs(
    x="PCA 1 (17.1%)",
    y="PCA 2 (5.9%)",
    title="PAO1 samples mapped to PAO1 reference",
)
fig1 += pn.scales.xlim(-50, 250)
fig1 += pn.scales.ylim(-60, 100)
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

print(fig1)

fig1.save("pao1_compendium_pao1_ref_pca.svg", format="svg", dpi=300)

# +
# Plot gene expression in PAO1 reference

# Plot pao1 and pa14 compendia separately to better visualize the overlap
fig2 = pn.ggplot(
    pao1_pca_encoded_label[pao1_pca_encoded_label["compendium"] == "PA14"],
    pn.aes(x="1", y="2"),
)
fig2 += pn.geom_point(color="#895881", alpha=0.3, size=4, stroke=0.8)
fig2 += pn.labs(
    x="PCA 1 (17.1%)",
    y="PCA 2 (5.9%)",
    title="PA14 samples mapped to PAO1 reference",
)
fig2 += pn.scales.xlim(-50, 250)
fig2 += pn.scales.ylim(-60, 100)
fig2 += pn.theme_bw()
fig2 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)

print(fig2)

fig2.save("pa14_compendium_pao1_ref_pca.svg", format="svg", dpi=300)

# +
# Plot gene expression in PA14 reference
fig3 = pn.ggplot(
    pa14_pca_encoded_label[pa14_pca_encoded_label["compendium"] == "PAO1"],
    pn.aes(x="1", y="2"),
)
fig3 += pn.geom_point(color="#C6A9B5", alpha=0.3, size=4, stroke=0.8)
fig3 += pn.labs(
    x="PCA 1 (16.2%)",
    y="PCA 2 (5.7%)",
    title="PAO1 samples mapped to PA14 reference",
)
fig3 += pn.scales.xlim(-50, 250)
fig3 += pn.scales.ylim(-60, 150)
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

print(fig3)

fig3.save("pao1_compendium_pa14_ref_pca.svg", format="svg", dpi=300)

# +
# Plot gene expression in PA14 reference
fig4 = pn.ggplot(
    pa14_pca_encoded_label[pa14_pca_encoded_label["compendium"] == "PA14"],
    pn.aes(x="1", y="2"),
)
fig4 += pn.geom_point(color="#895881", alpha=0.3, size=4, stroke=0.8)
fig4 += pn.labs(
    x="PCA 1 (16.2%)",
    y="PCA 2 (5.7%)",
    title="PA14 samples mapped to PA14 reference",
)
fig4 += pn.scales.xlim(-50, 250)
fig4 += pn.scales.ylim(-60, 150)
fig4 += pn.theme_bw()
fig4 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=15),
)

print(fig4)

fig4.save("pa14_compendium_pa14_ref_pca.svg", format="svg", dpi=300)
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

# ### Compare centroids

# Get PAO1 and PA14 sample ids
pao1_sample_ids = strain_labels.query("compendium=='PAO1'").index
pa14_sample_ids = strain_labels.query("compendium=='PA14'").index

# +
# Get separate embedding for PAO1, PA14 samples
pca_pao1 = PCA(n_components=200)
model_pca_pao1 = pca_pao1.fit(normalized_pao1_expression_df)

# Using PCA model generated using data aligned to PAO1 reference
encoded_pao1_pca_pao1_only = model_pca_pao1.transform(
    normalized_pao1_expression_df.loc[pao1_sample_ids]
)
encoded_pao1_pca_pa14_only = model_pca_pao1.transform(
    normalized_pao1_expression_df.loc[pa14_sample_ids]
)

# +
# Get separate embedding for PAO1, PA14 samples
pca_pa14 = PCA(n_components=200)
model_pca_pa14 = pca_pa14.fit(normalized_pa14_expression_df)

# Using PCA model generated using data aligned to PA14 reference
encoded_pa14_pca_pao1_only = model_pca_pa14.transform(
    normalized_pa14_expression_df.loc[pao1_sample_ids]
)
encoded_pa14_pca_pa14_only = model_pca_pa14.transform(
    normalized_pa14_expression_df.loc[pa14_sample_ids]
)

# +
# Get centroid using PAO1 reference
pao1_centroid_pao1 = np.mean(encoded_pao1_pca_pao1_only, axis=0)
pa14_centroid_pao1 = np.mean(encoded_pao1_pca_pa14_only, axis=0)

# We have a 200 dimensional array that is our centroid
print(pao1_centroid_pao1.shape)
print(pao1_centroid_pao1, pa14_centroid_pao1)

# +
# Get centroid using PA14 reference
pao1_centroid_pa14 = np.mean(encoded_pa14_pca_pao1_only, axis=0)
pa14_centroid_pa14 = np.mean(encoded_pa14_pca_pa14_only, axis=0)

print(pao1_centroid_pa14.shape)
print(pao1_centroid_pa14, pa14_centroid_pa14)

# +
both_centroid_pao1 = np.vstack((pao1_centroid_pao1, pa14_centroid_pao1))
both_centroid_pa14 = np.vstack((pao1_centroid_pa14, pa14_centroid_pa14))

print(both_centroid_pao1.shape)
print(both_centroid_pa14.shape)

# +
# Calculate distance between centroids
mean_dist_pao1 = pdist(both_centroid_pao1)
mean_dist_pa14 = pdist(both_centroid_pa14)


print("mean using PAO1 reference:", mean_dist_pao1)
print("mean using PA14 reference:", mean_dist_pa14)
# -

# ### Compare variance
# Here we are summing the variance along each PCs to represent the overall spread of the data (the spread along each direction)

# +
# Get variance of all samples using PAO1 reference
all_sample_ids = pao1_sample_ids.append(pa14_sample_ids)

encoded_pao1_pca = model_pca_pao1.transform(
    normalized_pao1_expression_df.loc[all_sample_ids]
)

encoded_pa14_pca = model_pca_pa14.transform(
    normalized_pa14_expression_df.loc[all_sample_ids]
)

print(encoded_pao1_pca.shape)
print(encoded_pa14_pca.shape)

# +
total_pao1_var = np.var(encoded_pao1_pca)
total_pa14_var = np.var(encoded_pa14_pca)

print("total PAO1 var: ", total_pao1_var)
print("total PA14 var: ", total_pa14_var)

# +
# Variance of PAO1, PA14 compendium (using PAO1 reference)
print(encoded_pao1_pca_pao1_only.shape)
print(encoded_pao1_pca_pa14_only.shape)

print(np.var(encoded_pao1_pca_pao1_only))
print(np.var(encoded_pao1_pca_pa14_only))

# +
# Variance of PAO1, PA14 compendium (using PA14 reference)
print(encoded_pa14_pca_pao1_only.shape)
print(encoded_pa14_pca_pa14_only.shape)

print(np.var(encoded_pa14_pca_pao1_only))
print(np.var(encoded_pa14_pca_pa14_only))
