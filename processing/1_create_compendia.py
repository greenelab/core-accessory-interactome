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
# A sample will be PA14 if:
# 1. PA14 mapping rate >= 30%
# 2. PA14-PAO1 mapping rate > 0%

import os
import pandas as pd
import seaborn as sns
from core_acc_modules import paths

# Params
mapping_threshold = 30
mapping_threshold_pa14 = 10
diff_mapping_threshold = 0

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

# ## TO DO: check why duplicates are appearing

# +
# Aggregate boolean labels into a single strain label
aggregated_label = []
for exp_id in list(sample_to_strain_table_full.index):
    if sample_to_strain_table_full.loc[exp_id, "PAO1"].all() == True:
        aggregated_label.append("PAO1")
    elif sample_to_strain_table_full.loc[exp_id, "PA14"].all() == True:
        aggregated_label.append("PA14")
    elif sample_to_strain_table_full.loc[exp_id, "PAK"].all() == True:
        aggregated_label.append("PAK")
    elif sample_to_strain_table_full.loc[exp_id, "ClinicalIsolate"].all() == True:
        aggregated_label.append("Clinical Isolate")
    else:
        aggregated_label.append("NA")

sample_to_strain_table_full["Strain type"] = aggregated_label

sample_to_strain_table = sample_to_strain_table_full["Strain type"].to_frame()
sample_to_strain_table = sample_to_strain_table.loc[
    sample_to_strain_table.index.drop_duplicates()
]

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
    pao1_logs.query("mapping_rate>=@mapping_threshold&diff_mapping_rate>0").index
)
high_pa14_mapping_ids = list(
    pa14_logs.query("mapping_rate>=@mapping_threshold&diff_mapping_rate>0").index
)

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

# +
# pao1_expression_label[pao1_expression_label.index.duplicated(keep=False)]

# +
# sample_to_strain_table[sample_to_strain_table.index.duplicated(keep=False)]
# -

print(pa14_expression_label.shape)
pa14_expression_label.head()

# ## Quick comparison
#
# Quick check comparing our binned labels compared with SRA annotations

pao1_expression_label["Strain type"].value_counts()

pa14_expression_label["Strain type"].value_counts()

# +
# Save compendia with label
pao1_expression_label.to_csv(paths.PAO1_COMPENDIUM_LABEL)
pa14_expression_label.to_csv(paths.PA14_COMPENDIUM_LABEL)

# Save compendia without label
pao1_expression_binned.to_csv(paths.PAO1_COMPENDIUM)
pa14_expression_binned.to_csv(paths.PA14_COMPENDIUM)

# Save metadata table
sample_to_strain_table.to_csv(paths.SAMPLE_TO_STRAIN_PROCESSED)
