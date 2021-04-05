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
# This notebook is using the thresholds from the [previous notebook](0_decide_thresholds.ipynb) to bin samples into PAO1 or PA14 compendia.
#
# A sample will be PAO1 if:
# 1. PAO1 mapping rate >= 30%
# 2. PAO1-PA14 mapping rate > 0%
#
# Note: if the difference in mapping rate is 0 then the same maps equally well to a PAO1 and PA14 reference. Here we looking for samples that preferentially map to PAO1.
#
# A sample will be PA14 if:
# 1. PA14 mapping rate >= 30%
# 2. PA14-PAO1 mapping rate > 0%

# %load_ext autoreload
# %autoreload 2
import os
import pandas as pd
import seaborn as sns
from textwrap import fill
import matplotlib.pyplot as plt
from core_acc_modules import paths

# Params
mapping_threshold = 30
diff_mapping_threshold = 2
diff_mapping_threshold_min = 0
diff_mapping_threshold_max = 2

# ## Load data

# +
# Log files
pao1_logs_filename = paths.PAO1_LOGS
pa14_logs_filename = paths.PA14_LOGS

# Expression data files
pao1_expression_filename = paths.PAO1_GE
pa14_expression_filename = paths.PA14_GE

# File containing table to map sample id to strain name
sample_to_strain_filename = paths.SAMPLE_TO_STRAIN

# +
# Load log files
pao1_logs = pd.read_csv(pao1_logs_filename, index_col=0, header=0)
pa14_logs = pd.read_csv(pa14_logs_filename, index_col=0, header=0)

# Load expression data
# Matrices will be sample x gene after taking the transpose
pao1_expression = pd.read_csv(pao1_expression_filename, index_col=0, header=0).T

pa14_expression = pd.read_csv(pa14_expression_filename, index_col=0, header=0).T

# Drop row with gene ensembl ids
pao1_expression.drop(["X"], inplace=True)
pa14_expression.drop(["X"], inplace=True)

# Load metadata
# Set index to experiment id, which is what we will use to map to expression data
sample_to_strain_table_full = pd.read_csv(sample_to_strain_filename, index_col=2)
# -

sample_to_strain_table_full.head()

pao1_logs.head()

pao1_expression.head()

# ## Format data
#
# Format index to only include experiment id. This will be used to map to expression data and labels

# +
# Format log indices so that values can be mapped to expression data
pao1_index_processed = pao1_logs.index.str.split("/").str[-1]
pa14_index_processed = pa14_logs.index.str.split("/").str[-1]

print(f"No. of samples processed using PAO1 reference: {pao1_logs.shape[0]}")
print(f"No. of samples processed using PA14 reference: {pa14_logs.shape[0]}")
pao1_logs.index = pao1_index_processed
pa14_logs.index = pa14_index_processed

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

pao1_logs.head()

pao1_expression.head()

# Save pre-binned expression data
pao1_expression.to_csv(paths.PAO1_PREBIN_COMPENDIUM, sep="\t")
pa14_expression.to_csv(paths.PA14_PREBIN_COMPENDIUM, sep="\t")

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

# ## Bin samples as PAO1 or PA14
#
# * Bin samples based on threshold from previous notebook
# * Check if there are any samples that have a high mapping to both PAO1 and PA14 (i.e. ambiguous mapping)

# Add column calculating the difference in mapping rates
pao1_logs["diff_mapping_rate"] = pao1_logs["mapping_rate"] - pa14_logs["mapping_rate"]
pa14_logs["diff_mapping_rate"] = pa14_logs["mapping_rate"] - pao1_logs["mapping_rate"]

# +
high_pao1_mapping_ids = list(
    pao1_logs.query(
        "mapping_rate>=@mapping_threshold&diff_mapping_rate>@diff_mapping_threshold"
    ).index
)
high_pa14_mapping_ids = list(
    pa14_logs.query(
        "mapping_rate>=@mapping_threshold&diff_mapping_rate>@diff_mapping_threshold"
    ).index
)

"""high_pao1_mapping_ids = list(
    pao1_logs.query("mapping_rate>=@mapping_threshold&diff_mapping_rate>@diff_mapping_threshold_min&diff_mapping_rate<@diff_mapping_threshold_max").index
)
high_pa14_mapping_ids = list(
    pa14_logs.query("mapping_rate>=@mapping_threshold&diff_mapping_rate>@diff_mapping_threshold_min&diff_mapping_rate<@diff_mapping_threshold_max").index
)"""

print(len(high_pao1_mapping_ids))
print(len(high_pa14_mapping_ids))

# +
# Check if any ids have high mapping rate for both PAO1 and PA14
high_pao1_pa14_mapping_ids = list(
    set(high_pao1_mapping_ids).intersection(high_pa14_mapping_ids)
)

print(len(high_pao1_pa14_mapping_ids))
# -

# **Some observations:**
# * Looks like there are not any samples that map to both PAO1 and PA14 using our criteria
# * The number of PA14 samples is much lower compared to PAO1. Does this mean that the mapping rates of PA14 samples mapped to PA14 reference lower?

# ## Create compendia
#
# Create PAO1 and PA14 compendia

# +
# Get expression data
# Note: reindexing needed here instead of .loc since samples from expression data
# were filtered out for low counts, but these samples still exist in log files
pao1_expression_binned = pao1_expression.reindex(high_pao1_mapping_ids)
pa14_expression_binned = pa14_expression.reindex(high_pa14_mapping_ids)

# Missing samples are dropped
pao1_expression_binned = pao1_expression_binned.dropna()
pa14_expression_binned = pa14_expression_binned.dropna()

# Drop ambiguously mapped samples
pao1_expression_binned = pao1_expression_binned.drop(high_pao1_pa14_mapping_ids)
pa14_expression_binned = pa14_expression_binned.drop(high_pao1_pa14_mapping_ids)
# -

print(pao1_expression_binned.shape)
print(pa14_expression_binned.shape)

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

# ## Quick comparison
#
# Quick check comparing our binned labels compared with SRA annotations

pao1_expression_label["Strain type"].value_counts()

pa14_expression_label["Strain type"].value_counts()

# ## Checks
#
# We noticed that the number of PA14 binned samples is much lower compared to the number of PAO1 samples. Let's look at the distribution of mapping rates for SRA annotated PAO1 and PA14 samples

# +
# Get SRA annotated sample ids
sra_pao1_sample_ids = list(
    sample_to_strain_table[sample_to_strain_table["Strain type"] == "PAO1"].index
)
sra_pa14_sample_ids = list(
    sample_to_strain_table[sample_to_strain_table["Strain type"] == "PA14"].index
)

print(len(sra_pao1_sample_ids))
print(len(sra_pa14_sample_ids))

# +
# Plot distribution of mapping rates to PAO1 and PA14

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8, 4))

# Distribution plot for core genes
sns.distplot(
    pao1_logs.loc[sra_pao1_sample_ids, "mapping_rate"],
    label="PAO1 (SRA annotated) mapping rate",
    color="red",
    kde=False,
    ax=axes[0],
)

sns.distplot(
    pa14_logs.loc[sra_pa14_sample_ids, "mapping_rate"],
    label="PA14 (SRA annotated) mapping rate",
    color="blue",
    kde=False,
    ax=axes[1],
)

plt.suptitle(
    fill("Distribution of mapping rates (SRA annotated)", width=40),
    x=0.5,
    y=1,
    fontsize=16,
)
axes[0].set_title(fill("PAO1 mapping rate", width=20))
axes[1].set_title(fill("PA14 mapping rate", width=20))
axes[0].set_xlabel("")
axes[1].set_xlabel("")
fig.text(0.5, 0.01, "Mapping rate", ha="center", fontsize=14)
fig.text(0.01, 0.5, "Count", ha="center", rotation=90, fontsize=14)
# -

# Looks like there are fewer PA14 samples with high PA14 mapping, which explains why we see such a reduced number of PA14 binned samples. We may need to used different thresholds for PAO1 and PA14.

# +
# Save compendia with label
pao1_expression_label.to_csv(paths.PAO1_COMPENDIUM_LABEL, sep="\t")
pa14_expression_label.to_csv(paths.PA14_COMPENDIUM_LABEL, sep="\t")

# Save compendia without label
pao1_expression_binned.to_csv(paths.PAO1_COMPENDIUM, sep="\t")
pa14_expression_binned.to_csv(paths.PA14_COMPENDIUM, sep="\t")

# Save metadata table
sample_to_strain_table.to_csv(paths.SAMPLE_TO_STRAIN_PROCESSED, sep="\t")
