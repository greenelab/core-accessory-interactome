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
#     display_name: Python [conda env:core_acc_env] *
#     language: python
#     name: conda-env-core_acc_env-py
# ---

# # Create PAO1 and PA14 compendia
#
# This notebook is using the observation from the [exploratory notebook](../explore_data/cluster_by_accessory_gene.ipynb) to bin samples into PAO1 or PA14 compendia.
#
# A sample is considered PAO1 if the median gene expression of PA14 accessory genes is 0 and PAO1 accessory genes in > 0.
# Similarlty, a sample is considered PA14 if the median gene expression of PA14 accessory genes is > 0 and PAO1 accessory genes in 0.

# %load_ext autoreload
# %autoreload 2
import os
import pandas as pd
import seaborn as sns
from textwrap import fill
import matplotlib.pyplot as plt
from core_acc_modules import paths, utils

# ## Load data

# +
# Expression data files
pao1_expression_filename = paths.PAO1_GE
pa14_expression_filename = paths.PA14_GE

# File containing table to map sample id to strain name
sample_to_strain_filename = paths.SAMPLE_TO_STRAIN

# +
# Load expression data
# Matrices will be sample x gene after taking the transpose
pao1_expression = pd.read_csv(pao1_expression_filename, index_col=0, header=0).T

pa14_expression = pd.read_csv(pa14_expression_filename, index_col=0, header=0).T

# Drop row with gene ensembl ids
pao1_expression.drop(["X"], inplace=True)
pa14_expression.drop(["X"], inplace=True)
# -

# Load metadata
# Set index to experiment id, which is what we will use to map to expression data
sample_to_strain_table_full = pd.read_csv(sample_to_strain_filename, index_col=2)

# ## Get core and accessory annotations

# +
pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT

core_acc_dict = utils.get_my_core_acc_genes(
    pao1_annot_filename, pa14_annot_filename, pao1_expression, pa14_expression
)
# -

pao1_acc = core_acc_dict["acc_pao1"]
pa14_acc = core_acc_dict["acc_pa14"]

# ## Format expression data
#
# Format index to only include experiment id. This will be used to map to expression data and SRA labels later

# +
# Format expression data indices so that values can be mapped to `sample_to_strain_table`
pao1_index_processed = pao1_expression.index.str.split(".").str[0]
pa14_index_processed = pa14_expression.index.str.split(".").str[0]

print(
    f"No. of samples processed using PAO1 reference after filtering: {pao1_expression.shape}"
)
print(
    f"No. of samples processed using PA14 reference after filtering: {pa14_expression.shape}"
)
pao1_expression.index = pao1_index_processed
pa14_expression.index = pa14_index_processed
# -

pao1_expression.head()

pa14_expression.head()

# Save pre-binned expression data
pao1_expression.to_csv(paths.PAO1_PREBIN_COMPENDIUM, sep="\t")
pa14_expression.to_csv(paths.PA14_PREBIN_COMPENDIUM, sep="\t")

# ## Bin samples as PAO1 or PA14

# +
# Create accessory df
# accessory gene ids | median accessory expression | strain label

# PAO1
pao1_acc_expression = pao1_expression[pao1_acc]
pao1_acc_expression["median_acc_expression"] = pao1_acc_expression.median(axis=1)

# PA14
pa14_acc_expression = pa14_expression[pa14_acc]
pa14_acc_expression["median_acc_expression"] = pa14_acc_expression.median(axis=1)

pao1_acc_expression.head()

# +
# Merge PAO1 and PA14 accessory dataframes
pao1_pa14_acc_expression = pao1_acc_expression.merge(
    pa14_acc_expression,
    left_index=True,
    right_index=True,
    suffixes=["_pao1", "_pa14"],
)

pao1_pa14_acc_expression.head()
# -

# Find PAO1 samples
pao1_binned_ids = list(
    pao1_pa14_acc_expression.query(
        "median_acc_expression_pao1>10 & median_acc_expression_pa14==0"
    ).index
)

# Find PA14 samples
pa14_binned_ids = list(
    pao1_pa14_acc_expression.query(
        "median_acc_expression_pao1==0 & median_acc_expression_pa14>10"
    ).index
)

# +
# Check that there are no samples that are binned as both PAO1 and PA14
shared_pao1_pa14_binned_ids = list(set(pao1_binned_ids).intersection(pa14_binned_ids))

assert len(shared_pao1_pa14_binned_ids) == 0
# -

# ## Format SRA annotations

# +
# Since experiments have multiple runs there are duplicated experiment ids in the index
# We will need to remove these so that the count calculations are accurate
sample_to_strain_table_full_processed = sample_to_strain_table_full[
    ~sample_to_strain_table_full.index.duplicated(keep="first")
]

assert (
    len(sample_to_strain_table_full.index.unique())
    == sample_to_strain_table_full_processed.shape[0]
)

# +
# Aggregate boolean labels into a single strain label
aggregated_label = []
for exp_id in list(sample_to_strain_table_full_processed.index):
    if sample_to_strain_table_full_processed.loc[exp_id, "PAO1"].all() == True:
        aggregated_label.append("PAO1")
    elif sample_to_strain_table_full_processed.loc[exp_id, "PA14"].all() == True:
        aggregated_label.append("PA14")
    elif sample_to_strain_table_full_processed.loc[exp_id, "PAK"].all() == True:
        aggregated_label.append("PAK")
    elif (
        sample_to_strain_table_full_processed.loc[exp_id, "ClinicalIsolate"].all()
        == True
    ):
        aggregated_label.append("Clinical Isolate")
    else:
        aggregated_label.append("NA")

sample_to_strain_table_full_processed["Strain type"] = aggregated_label

sample_to_strain_table = sample_to_strain_table_full_processed["Strain type"].to_frame()

sample_to_strain_table.head()
# -

# ## Create compendia
#
# Create PAO1 and PA14 compendia

# +
# Get expression data
# Note: reindexing needed here instead of .loc since samples from expression data
# were filtered out for low counts, but these samples still exist in log files
pao1_expression_binned = pao1_expression.loc[pao1_binned_ids]
pa14_expression_binned = pa14_expression.loc[pa14_binned_ids]

"""# Missing samples are dropped
pao1_expression_binned = pao1_expression_binned.dropna()
pa14_expression_binned = pa14_expression_binned.dropna()

# Drop ambiguously mapped samples
pao1_expression_binned = pao1_expression_binned.drop(high_pao1_pa14_mapping_ids)
pa14_expression_binned = pa14_expression_binned.drop(high_pao1_pa14_mapping_ids)"""
# -

assert len(pao1_binned_ids) == pao1_expression_binned.shape[0]
assert len(pa14_binned_ids) == pa14_expression_binned.shape[0]

# Label samples with SRA annotations
# pao1_expression_label = pao1_expression_binned.join(
#    sample_to_strain_table, how='left')
pao1_expression_label = pao1_expression_binned.merge(
    sample_to_strain_table, left_index=True, right_index=True
)
pa14_expression_label = pa14_expression_binned.merge(
    sample_to_strain_table, left_index=True, right_index=True
)
print(pao1_expression_label.shape)
pao1_expression_label.head()

print(pa14_expression_label.shape)
pa14_expression_label.head()

assert pao1_expression_binned.shape[0] == pao1_expression_label.shape[0]
assert pa14_expression_binned.shape[0] == pa14_expression_label.shape[0]

sample_to_strain_table["Strain type"].value_counts()

# Looks like our binned compendium sizes is fairly close in number to what SRA annotates

# ## Quick comparison
#
# Quick check comparing our binned labels compared with SRA annotations

pao1_expression_label["Strain type"].value_counts()

# **Manually check that these PA14 are mislabeled**
# * Clinical ones can be removed by increasing threshold

pa14_expression_label["Strain type"].value_counts()

# ## Check
#
# Manually look up the samples we binned as PAO1 but SRA labeled as PA14. Are these cases of samples being mislabeled?

pao1_expression_label[pao1_expression_label["Strain type"] == "PA14"]

# Note: These are the 7 PA14 labeled samples using threshold of 0
#
# Most samples appear to be mislabeled:
# * SRX5099522: https://www.ncbi.nlm.nih.gov/sra/?term=SRX5099522
# * SRX5099523: https://www.ncbi.nlm.nih.gov/sra/?term=SRX5099523
# * SRX5099524: https://www.ncbi.nlm.nih.gov/sra/?term=SRX5099524
# * SRX5290921: https://www.ncbi.nlm.nih.gov/sra/?term=SRX5290921
# * SRX5290922: https://www.ncbi.nlm.nih.gov/sra/?term=SRX5290922
#
# Two samples appear to be PA14 samples treated with antimicrobial manuka honey.
# * SRX7423386: https://www.ncbi.nlm.nih.gov/sra/?term=SRX7423386
# * SRX7423388: https://www.ncbi.nlm.nih.gov/sra/?term=SRX7423388

pa14_label_pao1_binned_ids = list(
    pao1_expression_label[pao1_expression_label["Strain type"] == "PA14"].index
)
pao1_pa14_acc_expression.loc[
    pa14_label_pao1_binned_ids,
    ["median_acc_expression_pao1", "median_acc_expression_pa14"],
]

# +
# Save compendia with SRA label
pao1_expression_label.to_csv(paths.PAO1_COMPENDIUM_LABEL, sep="\t")
pa14_expression_label.to_csv(paths.PA14_COMPENDIUM_LABEL, sep="\t")

# Save compendia without SRA label
pao1_expression_binned.to_csv(paths.PAO1_COMPENDIUM, sep="\t")
pa14_expression_binned.to_csv(paths.PA14_COMPENDIUM, sep="\t")

# Save processed metadata table
sample_to_strain_table.to_csv(paths.SAMPLE_TO_STRAIN_PROCESSED, sep="\t")
