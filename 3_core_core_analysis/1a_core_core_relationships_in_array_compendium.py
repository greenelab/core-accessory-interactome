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

# # Core-core relationships in array compendium
#
# This notebook performs the same stability analysis using the *P. aeruginosa* array compendium that is described in [Tan et al.](https://journals.asm.org/doi/10.1128/msystems.00025-15?permanently=true) and found in the associated repository [here](https://github.com/greenelab/adage/blob/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl). This notebook then compares the most stable genes identified using the array compendium compared to the RNA-seq compendium to validate our results are robust across platform (positive control).

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import scipy.stats
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as ss
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from scripts import utils, paths

random.seed(1)
# -

# Params
most_percent = 0.05
least_percent = 0.05

# +
# Output filenames
array_similarity_dist_filename = "array_similarity_scores_dist_spell.svg"

# Files containing genes with highest and lowest transcriptional similarity scores high and low
array_similarity_scores_filename = "array_similarity_scores_spell.tsv"

# +
# Import correlation matrix
array_expression_filename = paths.ARRAY_COMPENDIUM_GE
array_metadata_filename = paths.ARRAY_COMPENDIUM_METADATA

array_compendium = pd.read_csv(
    array_expression_filename, sep="\t", index_col=0, header=0
).T
array_metadata = pd.read_csv(array_metadata_filename, sep="\t", index_col=0, header=0)
# -

print(array_compendium.shape)
array_compendium.head()

print(array_metadata.shape)
array_metadata.head()

# ## Select only core genes

# +
# Read in RNA-seq transcriptional statistics
pao1_similarity_scores_filename = "pao1_core_similarity_associations_final_spell.tsv"

pao1_rnaseq_similarity_scores = pd.read_csv(
    pao1_similarity_scores_filename, sep="\t", header=0, index_col=0
)
# -

# Get only core gene ids
rnaseq_core_gene_ids = list(pao1_rnaseq_similarity_scores.index)
array_gene_ids = array_compendium.columns
array_core_gene_ids = set(rnaseq_core_gene_ids).intersection(array_gene_ids)

print(len(rnaseq_core_gene_ids))
print(len(array_gene_ids))
print(len(array_core_gene_ids))

array_compendium_core = array_compendium[array_core_gene_ids]

print(array_compendium_core.shape)
array_compendium_core.head()

# ## Make PAO1, PA14 array compendia

# Set "ml_data_source" (which corresponds to the sample ids in our expression matrix) as the index
array_metadata.set_index("ml_data_source", inplace=True)

# Select and separate between samples that are using PAO1 strain and those using PA14 strain
pao1_sample_ids = array_metadata.query("strain=='PAO1'").index.dropna()
pa14_sample_ids = array_metadata.query("strain=='PA14'").index.dropna()

print(len(pao1_sample_ids))
print(len(pa14_sample_ids))

# +
# Make PAO1, PA14 array compendia
pao1_sample_ids_shared = set(array_compendium_core.index).intersection(pao1_sample_ids)
pa14_sample_ids_shared = set(array_compendium_core.index).intersection(pa14_sample_ids)

pao1_array_compendium = array_compendium_core.loc[pao1_sample_ids_shared]
pa14_array_compendium = array_compendium_core.loc[pa14_sample_ids_shared]
# -

print(pao1_array_compendium.shape)
pao1_array_compendium.head()

print(pa14_array_compendium.shape)
pa14_array_compendium.head()

# ## Calculate correlation matrix
#
# Here we're following the same processing we performed for the RNA-seq data, to log10 transform the data and then apply SPELL.

num_SVs = 100

# Transpose compendia to be gene x sample
# Here we're interested in how genes cluster
pao1_array_compendium_T = pao1_array_compendium.T
pa14_array_compendium_T = pa14_array_compendium.T

# log transform data
pao1_array_compendium_log10 = np.log10(1 + pao1_array_compendium_T)
pa14_array_compendium_log10 = np.log10(1 + pa14_array_compendium_T)

# Apply SVD
pao1_U, pao1_s, pao1_Vh = np.linalg.svd(
    pao1_array_compendium_log10, full_matrices=False
)
pa14_U, pa14_s, pa14_Vh = np.linalg.svd(
    pa14_array_compendium_log10, full_matrices=False
)

print(pao1_array_compendium_T.shape)
print(pao1_U.shape, pao1_s.shape, pao1_Vh.shape)

print(pa14_array_compendium_T.shape)
print(pa14_U.shape, pa14_s.shape, pa14_Vh.shape)

# Convert ndarray to df to use corr()
pao1_U_df = pd.DataFrame(data=pao1_U, index=pao1_array_compendium_T.index)
pa14_U_df = pd.DataFrame(data=pa14_U, index=pa14_array_compendium_T.index)

# Correlation of U
# Since `corr()` computes pairwise correlation of columns we need to invert U
pao1_corr_log_spell = pao1_U_df.iloc[:, :num_SVs].T.corr()
pa14_corr_log_spell = pa14_U_df.iloc[:, :num_SVs].T.corr()

# Visually inspect that our log-SPELL processed data removes the dominant correlation signal that is present. Since the visualization can take a few minutes to plot, we've commented it out moving forward. Can refer to the saved plots

"""%%time
h1a = sns.clustermap(pao1_corr_log_spell, cmap="BrBG", center=0, figsize=(20, 20))
h1a.fig.suptitle(
    f"log transform + SPELL corrected using {num_SVs} vectors (PAO1)",
    y=1.05,
    fontsize=24,
)

# Save
pao1_log_spell_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_array_log_spell_clustermap.png"
)
h1a.savefig(pao1_log_spell_filename, dpi=300)"""

"""%%time
h2a = sns.clustermap(pa14_corr_log_spell, cmap="BrBG", center=0, figsize=(20, 20))
h2a.fig.suptitle(
    f"log transformed + SPELL corrected using {num_SVs} vectors (PA14)",
    y=1.05,
    fontsize=24,
)

# Save
pa14_log_spell_filename = os.path.join(
    paths.LOCAL_DATA_DIR, "pa14_array_log_spell_clustermap.png"
)
h2a.savefig(pa14_log_spell_filename, dpi=300)"""

# ## Calculate transcriptional stability
#
# All genes in this array compendium are core genes, since they were selected based on their hybridization so we don't need to map between PAO1 and PA14 ids types like we did for the RNA-seq analysis.

# +
rows = []
for gene_id in array_compendium_core.columns:
    # Make sure that the genes in the correlation profile are in the same order
    # in the PAO1 correlation matrix and the PA14 correlation matrix
    # Otherwise reorder
    if (pao1_corr_log_spell.index.equals(pa14_corr_log_spell.index)) & (
        pao1_corr_log_spell.columns.equals(pa14_corr_log_spell.columns)
    ):
        pass
    else:
        print("ordering is not the same, going to reorder...")
        pa14_corr_log_spell = pa14_corr_log_spell.loc[
            pao1_corr_log_spell.index, pao1_corr_log_spell.columns
        ]
    # Get correlation score
    # Make dataframe  with PAO1 id, PA14 homolog, correlation score
    corr_coef, pval = scipy.stats.pearsonr(
        pao1_corr_log_spell[gene_id], pa14_corr_log_spell[gene_id]
    )

    rows.append(
        {
            "PAO1 id": gene_id,
            "Transcriptional similarity across strains": corr_coef,
            "P-value": pval,
        }
    )

corr_summary_df = pd.DataFrame(rows)
# -

corr_summary_df.set_index("PAO1 id", inplace=True)

print(corr_summary_df.shape)
corr_summary_df.head()

# ## Plots

# Add label for most and least stable genes
corr_summary_df["label"] = ""

# Add label for most and least stable genes based on top X%
array_most_stable = corr_summary_df.sort_values(
    "Transcriptional similarity across strains", ascending=False
).head(round(most_percent * len(corr_summary_df)))
array_least_stable = corr_summary_df.sort_values(
    "Transcriptional similarity across strains", ascending=False
).tail(round(least_percent * len(corr_summary_df)))

array_most_threshold = array_most_stable.iloc[-1][
    "Transcriptional similarity across strains"
]
array_least_threshold = array_least_stable.iloc[0][
    "Transcriptional similarity across strains"
]
print(array_least_threshold, array_most_threshold)

corr_summary_df.loc[array_most_stable.index, "label"] = "most stable"
corr_summary_df.loc[array_least_stable.index, "label"] = "least stable"

# +
# Plot distribution of correlation scores
# This scores indicate how transcriptionally similar genes are across PAO1 and PA14 strains
fig_array = sns.displot(
    data=corr_summary_df,
    x="Transcriptional similarity across strains",
    hue="label",
    hue_order=["least stable", "most stable", ""],
    # label=["", "least stable", "most stable"],
    palette={"": "lightgrey", "least stable": "#a6aed0ff", "most stable": "#4e1c80"},
    legend=True,
    alpha=0.8,
    bins=np.linspace(0, 1, 50),
)
fig_array._legend.remove()

old_legend = fig_array._legend
handles = old_legend.legendHandles

legend = plt.legend(
    handles=[handles[0], handles[1]],
    labels=[
        fig_array._legend.texts[0].get_text(),
        fig_array._legend.texts[1].get_text(),
    ],
    bbox_to_anchor=(1.05, 0.6),
    loc="upper left",
    borderaxespad=0,
    fontsize=12,
)

plt.title("Stability of core genes across strain types", fontsize=14, y=1.1)
plt.xlabel("Transcriptional stability", fontsize=12)
plt.ylabel("Count", fontsize=12)
# -

# ## Compare most stable core genes
#
# We want to compare the most stable core genes obtained using the P. aeruginosa RNA-seq compendium vs the array compendium to validate our findings are robust.

pao1_rnaseq_similarity_scores_subset = pao1_rnaseq_similarity_scores[
    ["Transcriptional similarity across strains", "P-value", "Name", "label"]
]

all_similarity_scores = pao1_rnaseq_similarity_scores_subset.merge(
    corr_summary_df, left_index=True, right_index=True, suffixes=["_rnaseq", "_array"]
)

print(all_similarity_scores.shape)
all_similarity_scores.head()

# +
# Get most and least stable core genes
rnaseq_most_stable_genes = list(
    all_similarity_scores[all_similarity_scores["label_rnaseq"] == "most stable"].index
)
rnaseq_least_stable_genes = list(
    all_similarity_scores[all_similarity_scores["label_rnaseq"] == "least stable"].index
)

array_most_stable_genes = list(
    all_similarity_scores[all_similarity_scores["label_array"] == "most stable"].index
)
array_least_stable_genes = list(
    all_similarity_scores[all_similarity_scores["label_array"] == "least stable"].index
)

# +
# Compare
most_stable_venn = venn2(
    [set(rnaseq_most_stable_genes), set(array_most_stable_genes)],
    set_labels=("RNA-seq compendium", "Array compendium"),
)

most_stable_venn.get_patch_by_id("11").set_color("purple")
most_stable_venn.get_patch_by_id("11").set_edgecolor("none")
most_stable_venn.get_patch_by_id("11").set_alpha(0.3)
most_stable_venn.get_patch_by_id("01").set_color("blue")
most_stable_venn.get_patch_by_id("01").set_edgecolor("none")
most_stable_venn.get_patch_by_id("01").set_alpha(0.3)

plt.title("Most stable core genes", fontsize=16, fontname="Verdana")
for text in most_stable_venn.set_labels:
    text.set_fontsize(14)
    text.set_fontname("Verdana")

for text in most_stable_venn.subset_labels:
    text.set_fontsize(12)
    text.set_fontname("Verdana")

# Save figure
plt.savefig(
    "most_stable_array_vs_rnaseq_venn.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# +
least_stable_venn = venn2(
    [set(rnaseq_least_stable_genes), set(array_least_stable_genes)],
    set_labels=("RNA-seq compendium", "Array compendium"),
)

least_stable_venn.get_patch_by_id("11").set_color("purple")
least_stable_venn.get_patch_by_id("11").set_edgecolor("none")
least_stable_venn.get_patch_by_id("11").set_alpha(0.3)
least_stable_venn.get_patch_by_id("01").set_color("blue")
least_stable_venn.get_patch_by_id("01").set_edgecolor("none")
least_stable_venn.get_patch_by_id("01").set_alpha(0.3)

plt.title("Least stable core genes", fontsize=16, fontname="Verdana")
for text in least_stable_venn.set_labels:
    text.set_fontsize(14)
    text.set_fontname("Verdana")

for text in least_stable_venn.subset_labels:
    text.set_fontsize(12)
    text.set_fontname("Verdana")


# Save figure
plt.savefig(
    "least_stable_array_vs_rnaseq_venn",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# ### Plot correlation between transcriptional similarity values

# +
# Calculate correlation
r, p = stats.pearsonr(
    all_similarity_scores["Transcriptional similarity across strains_rnaseq"],
    all_similarity_scores["Transcriptional similarity across strains_array"],
)

print(r, p)

# +
# Plot correlation
fig = sns.jointplot(
    data=all_similarity_scores,
    x="Transcriptional similarity across strains_rnaseq",
    y="Transcriptional similarity across strains_array",
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes")

fig.set_axis_labels(
    "Transcriptional stability (RNA-seq)",
    "Transcriptional stability (Array)",
    fontsize=14,
    fontname="Verdana",
)
fig.fig.suptitle(
    "Stability RNA-seq vs Array", fontsize=16, fontname="Verdana", y=0.9, x=0.45
)

fig.savefig(
    "transcriptional_similarity_array_vs_rnaseq.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# +
# Get correlation for only those most stable genes found
# Select only those genes found to be most stable in either rna-seq or array compendia
most_stable_similarity_scores = all_similarity_scores[
    (all_similarity_scores["label_rnaseq"] == "most stable")
    | (all_similarity_scores["label_array"] == "most stable")
]

print(most_stable_similarity_scores.shape)
most_stable_similarity_scores.head()

# +
# Calculate correlation
r, p = stats.pearsonr(
    most_stable_similarity_scores["Transcriptional similarity across strains_rnaseq"],
    most_stable_similarity_scores["Transcriptional similarity across strains_array"],
)

print(r, p)

# +
# Plot correlation
fig = sns.jointplot(
    data=most_stable_similarity_scores,
    x="Transcriptional similarity across strains_rnaseq",
    y="Transcriptional similarity across strains_array",
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes")

fig.set_axis_labels(
    "Transcriptional similarity (RNA-seq)",
    "Transcriptional similarity (Array)",
    fontsize=14,
    fontname="Verdana",
)
fig.fig.suptitle(
    "Stability RNA-seq vs Array (most stable only)",
    fontsize=16,
    fontname="Verdana",
    y=0.9,
    x=0.45,
)

fig.savefig(
    "transcriptional_similarity_array_vs_rnaseq_most_stable.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# Hypergeometric test:
#
# Given $N$ (`total_num_genes`) number of genes with $K$ (`num_rnaseq_stable`) most stable genes in RNA-seq compendium. Using the array compendium, we identify genes as being most stable. What is the probability that  of the genes identified in the array compendium are also most stable in the RNA-seq compendium? What is the probability of drawing $k$ (`num_concordant_stable_genes`) or more concordant genes?
#
# This was a way for us to quantify the correlation between RNA-seq and array most stable findings.

# Try calculating the enrichment of most stable core genes found from rna-seq in array
total_num_genes = all_similarity_scores.shape[0]
num_concordant_stable_genes = all_similarity_scores[
    (all_similarity_scores["label_rnaseq"] == "most stable")
    & (all_similarity_scores["label_array"] == "most stable")
].shape[0]
num_rnaseq_stable_genes = all_similarity_scores[
    all_similarity_scores["label_rnaseq"] == "most stable"
].shape[0]
num_array_stable_genes = all_similarity_scores[
    all_similarity_scores["label_array"] == "most stable"
].shape[0]

print(total_num_genes)
print(num_rnaseq_stable_genes)
print(num_array_stable_genes)
print(num_concordant_stable_genes)

p = ss.hypergeom.sf(
    num_concordant_stable_genes,
    total_num_genes,
    num_rnaseq_stable_genes,
    num_array_stable_genes,
)
print(p)

# ## Examine genes that differ

# +
most_rnaseq_only = set(rnaseq_most_stable_genes).difference(array_most_stable_genes)
most_array_only = set(array_most_stable_genes).difference(rnaseq_most_stable_genes)

least_rnaseq_only = set(rnaseq_least_stable_genes).difference(array_least_stable_genes)
least_array_only = set(array_least_stable_genes).difference(rnaseq_least_stable_genes)
# -

all_similarity_scores.loc[most_rnaseq_only]

all_similarity_scores.loc[most_array_only]

sns.displot(
    all_similarity_scores.loc[
        most_array_only, "Transcriptional similarity across strains_rnaseq"
    ]
)

all_similarity_scores.loc[least_rnaseq_only]

all_similarity_scores.loc[least_array_only]

# Save
fig_array.savefig(
    array_similarity_dist_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# **Takeaway:**
# * Venn diagram comparing the most stable core genes using the array vs RNA-seq compendium, where the most stable genes = top 5% of genes with the highest transcriptional similarity.
# * There is some consistency of most stable genes using the array and RNA-seq compendia
# * Looking at the RNA-seq transcriptional similarity of the array only genes (191 genes), most of these genes were just below the threshold for "most stable" (transcriptional similarity = 0.4 - 0.55) according to the distribution plot. However there are some genes that have a lower similarity score (0.2 - 0.3). Why are these stable in microarray but not RNA-seq?
#
# * I also plotted the results comparing the least stable genes across the compendium, but its less clear what this is telling us since these genes are those that have unstable profiles across strains.

# Save transcriptional similarity df
corr_summary_df.to_csv(array_similarity_scores_filename, sep="\t")