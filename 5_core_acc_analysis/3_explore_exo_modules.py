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

# # Explore exoS/exoU genes
#
# This notebook specifically explores the core genes related to the two exoS (PAO1) and exoU (PA14) accessory genes. Specifically examining the core genes that are highly co-expressed with both exoS and exoU versus those that are specific to one gene or the other.
#
# _P. aeruginosa_ uses a type III secretion system (T3SS) to promote development of severe disease, particularly in patients with impaired immune defenses. _P. aeruginosa_ uses a type III secretion system to inject toxic effector proteins into the cytoplasm of eukaryotic cells. ExoU, ExoS, and ExoT, three effector proteins secreted by this system. ExoU and ExoS are usually secreted by different strains.
#
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC529154/

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
from scripts import utils, paths, annotations

np.random.seed(1)
# -

# Get gene id for exoS and exoU
exoS_id = "PA3841"
exoU_id = "PA14_51530"

# ### Get relationship between exoS/U and all other genes

# +
# Read in correlation for all genes
pao1_all_corr_filename = paths.PAO1_CORR_LOG_SPELL
pa14_all_corr_filename = paths.PA14_CORR_LOG_SPELL

pao1_all_corr = pd.read_csv(pao1_all_corr_filename, sep="\t", index_col=0, header=0)
pa14_all_corr = pd.read_csv(pa14_all_corr_filename, sep="\t", index_col=0, header=0)
# -

# Get correlation between exoS/U and all other genes
exoS_all_corr = pao1_all_corr.loc[exoS_id].to_frame("corr to exoS")
exoU_all_corr = pa14_all_corr.loc[exoU_id].to_frame("corr to exoU")

print(exoS_all_corr.shape)
exoS_all_corr.head()

print(exoU_all_corr.shape)
exoU_all_corr.head()

# ### Add gene name

# +
# Import gene metadata
pao1_gene_annot_filename = paths.GENE_PAO1_ANNOT
pa14_gene_annot_filename = paths.GENE_PA14_ANNOT

pao1_gene_annot = pd.read_csv(pao1_gene_annot_filename, index_col=0, header=0)
pa14_gene_annot = pd.read_csv(pa14_gene_annot_filename, index_col=0, header=0)
# -

# Get df with gene ids as indices and gene names as a column
# Having the data in a df instead of a series will just allow me to do my merges that are in the notebook
pao1_gene_annot = pao1_gene_annot["Name"].to_frame("gene name")
pa14_gene_annot = pa14_gene_annot["Name"].to_frame("gene name")

# Add gene name annotation
exoS_all_corr = exoS_all_corr.merge(
    pao1_gene_annot, left_index=True, right_index=True, how="left"
)
exoU_all_corr = exoU_all_corr.merge(
    pa14_gene_annot, left_index=True, right_index=True, how="left"
)

print(exoS_all_corr.shape)
exoS_all_corr.head()

print(exoU_all_corr.shape)
exoU_all_corr.head()

# ### Select only core genes

# +
# Get core genes
pao1_expression_filename = paths.PAO1_COMPENDIUM
pa14_expression_filename = paths.PA14_COMPENDIUM

pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT
# -

# Make a dataframe with gene ids
pao1_gene_ids = pd.DataFrame(data=[], index=pao1_all_corr.index)
pa14_gene_ids = pd.DataFrame(data=[], index=pa14_all_corr.index)

(
    pao1_arr,
    pa14_arr,
    pao1_core,
    pao1_acc,
    pa14_core,
    pa14_acc,
) = annotations.map_core_acc_annot(
    pao1_gene_ids,
    pa14_gene_ids,
    pao1_expression_filename,
    pa14_expression_filename,
    pao1_annot_filename,
    pa14_annot_filename,
)

exoS_core_corr = exoS_all_corr.loc[pao1_core]
exoU_core_corr = exoU_all_corr.loc[pa14_core]

print(exoS_core_corr.shape)
print(exoU_core_corr.shape)

exoS_core_corr.head()

exoU_core_corr.head()

# ### Map and merge dataframes

gene_mapping_pa14 = utils.get_pao1_pa14_gene_map(pa14_annot_filename, "pa14")

pa14_gene_name_map = gene_mapping_pa14["PAO1_ID"].to_dict()

# Map PA14 gene ids to PAO1
exoU_core_corr = exoU_core_corr.rename(mapper=pa14_gene_name_map, axis=0)

print(exoU_core_corr.shape)
exoU_core_corr.head()

# +
# Merge dataframes to get core genes related to exoS and exoU in one dataframe
exo_core_corr = exoS_core_corr.merge(
    exoU_core_corr, left_index=True, right_index=True, how="inner"
)

print(exo_core_corr.shape)
exo_core_corr.head()
# -

# ### Plot

# Core genes highly co-expressed with both exoS and exoU
exo_core_both = exo_core_corr[
    (exo_core_corr["corr to exoS"] > 0.4) & (exo_core_corr["corr to exoU"] > 0.2)
]
exo_core_both_ids = exo_core_both.index
exo_core_both

# Core genes co-expressed with exoS
exoS_core_only = exo_core_corr[
    (exo_core_corr["corr to exoS"] > 0.2) & (exo_core_corr["corr to exoU"] < 0)
]
exoS_core_only_ids = exoS_core_only.index
exoS_core_only

# Add labels
exo_core_corr["label"] = ""
exo_core_corr.loc[exo_core_both_ids, "label"] = "both"
exo_core_corr.loc[exoS_core_only_ids, "label"] = "exoS only"

# +
plt.figure(figsize=[10, 8])
fig_exo_corr = sns.scatterplot(
    data=exo_core_corr,
    x="corr to exoS",
    y="corr to exoU",
    alpha=0.7,
    hue="label",
    palette={"": "darkgrey", "both": "#fd5e0c", "exoS only": "#f9da76"},
)

plt.ylabel(r"Correlation to $exoU$", fontsize=14)
plt.xlabel(R"Correlation to $exoS$", fontsize=14)
plt.title("Correlation of core genes with T3SS accessory genes", fontsize=16)
plt.legend(bbox_to_anchor=(1.05, 1), fontsize=14)
# -

sns.jointplot(data=exo_core_corr, x="corr to exoS", y="corr to exoU", kind="hex")

# +
# Save
exo_core_corr.to_csv("core_genes_related_to_exoSU.tsv", sep="\t")

fig_exo_corr.figure.savefig(
    "core_genes_correlated_with_exo.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# **Takeaway**
# * Core genes that are highly co-expressed with both exoS and exoU are related to the T3SS secretion machinery
# * Core genes highly co-expressed with exoS are TBD
