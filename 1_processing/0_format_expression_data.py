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

# # Format expression
#
# This notebook cleans and formats the data. There are two specific things that this notebook does:
# 1. We noticed that there are PAO1 gene ids in the features/column headers for PA14 gene expression matrix. This notebook removes those gene ids before we perform the rest of the downstream analyses.
# 2. This notebook also formats the data to be sample x gene matrices

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
from scripts import paths

# Raw (normalized counts) expression data files
pao1_expression_filename = paths.PAO1_GE_RAW
pa14_expression_filename = paths.PA14_GE_RAW

# Load expression data
# Matrices will be sample x gene after taking the transpose
pao1_expression = pd.read_csv(pao1_expression_filename, index_col=0, header=0).T
pa14_expression = pd.read_csv(pa14_expression_filename, index_col=0, header=0).T

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

print(pao1_expression.shape)
pao1_expression.head()

print(pa14_expression.shape)
pa14_expression.head()

# ## Clean
#
# Find any gene ids that are mismatched

for gene_id in pao1_expression.columns:
    if "PA14_" in gene_id:
        print(gene_id)

mismatched_gene_ids = []
for gene_id in pa14_expression.columns:
    if "PA14_" not in gene_id:
        print(gene_id)
        mismatched_gene_ids.append(gene_id)

# Drop this columns from the PA14 expression data
pa14_expression_new = pa14_expression.drop(columns=mismatched_gene_ids)

print(pa14_expression_new.shape)
pa14_expression_new.head()

# Verify that we've removed the gene ids
assert pa14_expression_new.shape[1] == pa14_expression.shape[1] - len(
    mismatched_gene_ids
)
assert pa14_expression_new.shape[0] == pa14_expression.shape[0]

for gene_id in mismatched_gene_ids:
    if gene_id in pa14_expression_new.columns:
        print(gene_id)

# ## Save

pao1_expression.to_csv(paths.PAO1_GE, sep="\t")
pa14_expression_new.to_csv(paths.PA14_GE, sep="\t")
