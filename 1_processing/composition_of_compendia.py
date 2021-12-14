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

# # Composition of compendia
#
# This notebook makes figures that illustrate the composition of the compendia. In particular this notebook takes metadata files (generated by Sam Neff from Dartmouth) and creates plots to show the different types of media used in experiments and what types of genetic malnipulations were used as well.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import plotnine as pn
import seaborn as sns
import matplotlib.pyplot as plt

# Import metadata files
pao1_metadata_filename = "PAO1TableVF1.csv"
pa14_metadata_filename = "PA14TableVF1.csv"

pao1_metadata = pd.read_csv(pao1_metadata_filename, header=0, index_col=0)
pa14_metadata = pd.read_csv(pa14_metadata_filename, header=0, index_col=0)

print(pao1_metadata.shape)
pao1_metadata.head(10)

print(pa14_metadata.shape)
pa14_metadata.head(10)

# +
# TO DO
# Clean up values?
# Only plot the top 10?
# Coloring palette
# -

# Format dataframe
# Only keep first row
pao1_metadata_first = pao1_metadata[~pao1_metadata.index.duplicated(keep="first")]
pa14_metadata_first = pa14_metadata[~pa14_metadata.index.duplicated(keep="first")]

# Concatenate
both_metadata_first = pd.concat([pao1_metadata_first, pa14_metadata_first])

# Rename to clean up legend
media_mapper = {
    "TSB": "Tryptic Soy Broth",
}
# gene_function_mapper = {}
both_metadata_first.replace({"Medium": media_mapper}, inplace=True)

# ## Plot media distribution in PAO1 and PA14
#
# Media will be at the study level

both_metadata_first_media = (
    both_metadata_first.groupby(["Strain", "Medium"])
    .size()
    .reset_index()
    .pivot(columns="Medium", index="Strain", values=0)
)

both_metadata_first_media

fig_media = both_metadata_first_media.plot(
    kind="bar", stacked=True, colormap="Set2", figsize=(12, 8)
)
plt.legend(bbox_to_anchor=(1.5, 1), loc="upper right", ncol=1)
plt.title("Media used in experiments", fontsize=16)
fig_media.set_xlabel("")
fig_media.set_ylabel("Count", fontsize=14)
fig_media.set_xticklabels(
    ["PA14 compendium", "PAO1 compendium"], rotation=0, fontsize=14
)

# ## Plot Gene function distribution in PAO1 and PA14
# Gene function will be at the study level as well since the gene will be the same

# Format
both_metadata_first_function = (
    both_metadata_first.groupby(["Strain", "Gene.Function"])
    .size()
    .reset_index()
    .pivot(columns="Gene.Function", index="Strain", values=0)
)

both_metadata_first_function

fig_function = both_metadata_first_function.plot(
    kind="bar", stacked=True, colormap="Set2", figsize=(12, 8)
)
plt.legend(bbox_to_anchor=(1.8, 1), loc="upper right", ncol=1)
plt.title("Gene function studied in experiments", fontsize=16)
fig_function.set_xlabel("")
fig_function.set_ylabel("Count", fontsize=14)
fig_function.set_xticklabels(
    ["PA14 compendium", "PAO1 compendium"], rotation=0, fontsize=14
)

# ## Plot pathways associated with perturbed gene
#
# There are 3 possible reasons that the value = Nan:
# 1. The sample was a control sample
# 2. The perturbed gene was not found to be associated to a pathway due to limitations in annotation information
# 3. The experimental design was not performing a genetic malnipulation

# +
# Format dataframe so that the only row with NaN corresponds to study where
# they did not perform a genetic malnipulation and so "Perturbed.Gene" = None
# -

# Concatenate
both_metadata_all = pd.concat([pao1_metadata, pa14_metadata])

both_metadata_kegg = (
    both_metadata_all.groupby(["Strain", "KEGG.Pathway"])
    .size()
    .reset_index()
    .pivot(columns="KEGG.Pathway", index="Strain", values=0)
)

both_metadata_kegg

fig_kegg = both_metadata_kegg.plot(
    kind="bar", stacked=True, colormap="Set2", figsize=(12, 8)
)
plt.legend(bbox_to_anchor=(1.45, 1), loc="upper right", ncol=1)
plt.title("KEGG pathway studied in experiments", fontsize=16)
fig_kegg.set_xlabel("")
fig_kegg.set_ylabel("Count", fontsize=14)
fig_kegg.set_xticklabels(
    ["PA14 compendium", "PAO1 compendium"], rotation=0, fontsize=14
)

# Save plots
fig_media.figure.savefig(
    "compendia_media.svg", dpi=300, format="svg", bbox_inches="tight"
)
fig_function.figure.savefig(
    "compendia_gene_function.svg", dpi=300, format="svg", bbox_inches="tight"
)
fig_kegg.figure.savefig(
    "compendia_kegg.svg", dpi=300, format="svg", bbox_inches="tight"
)
