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

# # Process array data
#
# This notebook will process array data to:
# 1. Include only those PAO1 samples, based on the metadata
# 2. Use only genes shared by both the array and RNA-seq compendia

# %load_ext autoreload
# %autoreload 2
import pandas as pd
from core_acc_modules import paths

# +
# Files
array_compendium_filename = paths.ARRAY_DATA
array_metadata_filename = paths.ARRAY_METADATA

pao1_rnaseq_compendium_filename = paths.PAO1_COMPENDIUM
# -

array_compendium = pd.read_csv(
    array_compendium_filename, sep="\t", index_col=0, header=0
).T
array_metadata = pd.read_csv(array_metadata_filename, sep="\t", index_col=0, header=0)

print(array_compendium.shape)
array_compendium.head()

print(array_metadata.shape)
array_metadata.head()

pao1_strain_values = [
    strain_name
    for strain_name in array_metadata["strain"].unique()
    if ("PAO1" in strain_name) & ("MPAO1" not in strain_name)
]

pao1_sample_ids = array_metadata[
    array_metadata["strain"].isin(pao1_strain_values)
].loc[:, "ml_data_source"]

# Drop any sample ids that are na
pao1_sample_ids.dropna(inplace=True)

# ## Select subset of samples
#
# Select only those sample ids that are using PAO1-like strains and have expression data available

pao1_array_compendium = array_compendium.loc[set(pao1_sample_ids.values).intersection(array_compendium.index)]

print(pao1_array_compendium.shape)
pao1_array_compendium.head()

# Drop samples without expression data available
pao1_array_compendium.dropna(inplace=True)

print(pao1_array_compendium.shape)
pao1_array_compendium.head()

# ## Use only shared genes
#
# Use only genes that are shared between the array and RNA-seq compendia so that we can compare modules

# Read in rnaseq compendium processing
pao1_rnaseq_compendium = pd.read_csv(
    pao1_rnaseq_compendium_filename, sep="\t", index_col=0, header=0
)
print(pao1_rnaseq_compendium.shape)
pao1_rnaseq_compendium.head()

# +
# Get shared genes
pao1_array_gene_ids = pao1_array_compendium.columns
pao1_rnaseq_gene_ids = pao1_rnaseq_compendium.columns

shared_gene_ids = list(set(pao1_array_gene_ids).intersection(set(pao1_rnaseq_gene_ids)))
print(len(shared_gene_ids))
# -

# Only include shared genes for both compendia
pao1_array_compendium_processed = pao1_array_compendium[shared_gene_ids]
pao1_rnaseq_compendium_processed = pao1_rnaseq_compendium[shared_gene_ids]

# Save to new paths
pao1_array_compendium_processed.to_csv(paths.ARRAY_COMPENDIUM_TO_COMPARE, sep="\t")
pao1_rnaseq_compendium_processed.to_csv(paths.RNASEQ_COMPENDIUM_TO_COMPARE, sep="\t")
