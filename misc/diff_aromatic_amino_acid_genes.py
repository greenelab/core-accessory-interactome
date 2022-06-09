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

# # Explore differences in aromatic amino acid genes
#
# This notebook compares the expression activity of genes related to the aromatic amino acid pathway between PAO1 and PA14 strains using both RNA-seq and array compendia. The goal is to determine if the differences in expression activity between PAO1 vs PA14 are consistent using both RNA-seq and array compendia.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from scripts import paths

# ## Load expression compendia

# +
# RNA-seq compendia files
rnaseq_pao1_compendium_filename = paths.PAO1_COMPENDIUM
rnaseq_pa14_compendium_filename = paths.PA14_COMPENDIUM

# Array compendia files
array_pao1_compendium_filename = paths.PAO1_COMPENDIUM_ARRAY
array_pa14_compendium_filename = paths.PA14_COMPENDIUM_ARRAY

# +
# Load expression data
rnaseq_pao1_compendium = pd.read_csv(
    rnaseq_pao1_compendium_filename, sep="\t", index_col=0, header=0
)
rnaseq_pa14_compendium = pd.read_csv(
    rnaseq_pa14_compendium_filename, sep="\t", index_col=0, header=0
)

array_pao1_compendium = pd.read_csv(
    array_pao1_compendium_filename, sep="\t", index_col=0, header=0
)
array_pa14_compendium = pd.read_csv(
    array_pa14_compendium_filename, sep="\t", index_col=0, header=0
)
# -

# ## Load KEGG pathways

# Load KEGG pathway data
pao1_pathway_filename = "../3_core_core_analysis/pao1_kegg_annot.tsv"
pa14_pathway_filename = "../3_core_core_analysis/pa14_kegg_annot.tsv"

pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", header=0, index_col=0)
pa14_pathways = pd.read_csv(pa14_pathway_filename, sep="\t", header=0, index_col=0)

print(pao1_pathways.shape)
pao1_pathways.head()

print(pa14_pathways.shape)
pa14_pathways.head()

# ## Get aromatic amino acid genes

# +
# Aromatic amino acid pathways
# Are there others to include?
# https://docs.google.com/spreadsheets/d/1pDehBGTyCN0hTrlvbBNymeI640Sl68XlUBwBU8aPZTA/edit#gid=2100140648
pao1_aa_pathways = [
    "path:pae00362 : Benzoate degradation",
    "path:pae00350 : Tyrosine metabolism",
    "path:pae00620 : Pyruvate metabolism",
    "path:pae01220 : Degradation of aromatic compounds",
]

pa14_aa_pathways = [
    "path:pau00362 : Benzoate degradation",
    "path:pau00350 : Tyrosine metabolism",
    "path:pau00620 : Pyruvate metabolism",
    "path:pau01220 : Degradation of aromatic compounds",
]


# -

# Select associated gene sets
def get_associated_genes(pathway_db, selected_pathways):
    aa_genes = set()
    for gene_set in pathway_db.loc[selected_pathways, "gene_ids"]:
        gene_set_processed = set(json.loads(gene_set.replace("'", '"')))
        aa_genes = aa_genes.union(gene_set_processed)

    return list(aa_genes)


pao1_aa_genes = get_associated_genes(pao1_pathways, pao1_aa_pathways)
print(pao1_aa_genes, len(pao1_aa_genes))

pa14_aa_genes = get_associated_genes(pa14_pathways, pa14_aa_pathways)
print(pa14_aa_genes, len(pa14_aa_genes))

# ## Get expression activity of aromatic amino acid genes

# +
# Select rna-seq genes that are related to AA pathways
rnaseq_pao1_aa_expression = rnaseq_pao1_compendium[pao1_aa_genes]
rnaseq_pa14_aa_expression = rnaseq_pa14_compendium[pa14_aa_genes]

rnaseq_pao1_aa_expression.head()

# +
# Select array genes that are related to AA pathways
shared_pao1_aa_genes = set(pao1_aa_genes).intersection(array_pao1_compendium.columns)
array_pao1_aa_expression = array_pao1_compendium[shared_pao1_aa_genes]
array_pa14_aa_expression = array_pa14_compendium[shared_pao1_aa_genes]

array_pao1_aa_expression.head()

# +
# Calculate median across samples
rnaseq_pao1_aa_expression_median = rnaseq_pao1_aa_expression.median()
rnaseq_pa14_aa_expression_median = rnaseq_pa14_aa_expression.median()

array_pao1_aa_expression_median = array_pao1_aa_expression.median()
array_pa14_aa_expression_median = array_pa14_aa_expression.median()
# -

# Calculate difference in median of medians
diff_rnaseq = abs(
    rnaseq_pao1_aa_expression_median.median()
    - rnaseq_pa14_aa_expression_median.median()
)
diff_rnaseq

diff_array = abs(
    array_pao1_aa_expression_median.median() - array_pa14_aa_expression_median.median()
)
diff_array

# ## Plot

# Plot box plot of rna-seq and array
rnaseq_pao1_aa_expression_median.plot(kind="density", label="PAO1", legend=True)
rnaseq_pa14_aa_expression_median.plot(kind="density", label="PA14", legend=True)
plt.title("Median expression of AA genes using RNA-seq")
plt.xlabel("Median expression")

array_pao1_aa_expression_median.plot(kind="density", label="PAO1", legend=True)
array_pa14_aa_expression_median.plot(kind="density", label="PA14", legend=True)
plt.title("Median expression of AA genes using array")
plt.xlabel("Median expression")
