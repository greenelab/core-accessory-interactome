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

# # Add annotations
#
# This notebook takes the dataframe with information about module composition and their labels and adds additional annotations including:
#
# 1. Which gene is contained within the modules (both gene id and gene name)
# 2. KEGG pathways that genes are found in
# 3. GO pathways genes are found in
# 4. Regulon/operon genes are found in

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from core_acc_modules import utils, paths

random.seed(1)
# -

# User param
method = "affinity"

# +
# Import module labels
pao1_module_label_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_gene_module_labels_{method}.tsv"
)
pa14_module_label_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_gene_module_labels_{method}.tsv"
)

pao1_module_labels = pd.read_csv(
    pao1_module_label_filename, sep="\t", index_col=0, header=0
)
pa14_module_labels = pd.read_csv(
    pa14_module_label_filename, sep="\t", index_col=0, header=0
)

# +
# Import gene memberships
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{method}.tsv"
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", index_col=0, header=0)
pa14_membership = pd.read_csv(pa14_membership_filename, sep="\t", index_col=0, header=0)

# +
# Import gene metadata
pao1_gene_annot_filename = paths.GENE_PAO1_ANNOT
pa14_gene_annot_filename = paths.GENE_PA14_ANNOT

pao1_gene_annot = pd.read_csv(pao1_gene_annot_filename, index_col=0, header=0)
pa14_gene_annot = pd.read_csv(pa14_gene_annot_filename, index_col=0, header=0)
# -

pao1_gene_annot = pao1_gene_annot["Name"].to_frame("gene name")
pa14_gene_annot = pa14_gene_annot["Name"].to_frame("gene name")

# ## Add module labels

# Add module labels
pao1_gene_module_labels = pao1_membership.merge(
    pao1_module_labels, left_on="module id", right_index=True
)
pa14_gene_module_labels = pa14_membership.merge(
    pa14_module_labels, left_on="module id", right_index=True
)

# ## Add gene names

# Add gene names
pao1_gene_module_labels = pao1_gene_module_labels.merge(
    pao1_gene_annot, left_index=True, right_index=True
)
pa14_gene_module_labels = pa14_gene_module_labels.merge(
    pa14_gene_annot, left_index=True, right_index=True
)

print(pao1_gene_module_labels.shape)
pao1_gene_module_labels.head()

print(pa14_gene_module_labels.shape)
pa14_gene_module_labels.head()

# ## Add core/accessory annotations

# +
# Read in expression data
pao1_expression_filename = paths.PAO1_COMPENDIUM
pa14_expression_filename = paths.PA14_COMPENDIUM

pao1_expression = pd.read_csv(pao1_expression_filename, sep="\t", index_col=0, header=0)
pa14_expression = pd.read_csv(pa14_expression_filename, sep="\t", index_col=0, header=0)
# -

core_acc_dict = utils.get_my_core_acc_genes(
    pao1_gene_annot_filename, pa14_gene_annot_filename, pao1_expression, pa14_expression
)

pao1_core = core_acc_dict["core_pao1"]
pa14_core = core_acc_dict["core_pa14"]
pao1_acc = core_acc_dict["acc_pao1"]
pa14_acc = core_acc_dict["acc_pa14"]

pao1_gene_module_labels.loc[pao1_core, "core/acc"] = "core"
pao1_gene_module_labels.loc[pao1_acc, "core/acc"] = "acc"

pa14_acc_shared = set(pa14_acc).intersection(pa14_gene_module_labels.index)
pa14_gene_module_labels.loc[pa14_core, "core/acc"] = "core"
pa14_gene_module_labels.loc[pa14_acc_shared, "core/acc"] = "acc"

pao1_gene_module_labels.head()

pa14_gene_module_labels.head()

# ## Add KEGG pathways
#
# For each pathway, what genes are contained in it

# +
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", index_col=0, header=None)
# -

pao1_pathways[2] = pao1_pathways[2].str.split(";").apply(set)
pao1_pathways.index = pao1_pathways.index.str.split(" - ").str[0]
pao1_pathways.head()

gene_to_pathways_df = pd.DataFrame(
    index=pao1_gene_module_labels.index, columns=list(pao1_pathways.index)
)

# %%time
for gene in gene_to_pathways_df.index:
    gene_to_pathways_df.loc[gene] = [
        gene in pao1_pathways.loc[pathway, 2] for pathway in pao1_pathways.index
    ]

# ## Add operon
#
# For each operon, what genes are contained in it
#
# NOTE: This code takes a while to run so for now its commented out

# +
pao1_operon_filename = "https://raw.githubusercontent.com/greenelab/core-accessory-interactome/6635c0e357c0172c2cebd0368648030e0ee4beaf/data/metadata/operons_format.csv"

pao1_operons = pd.read_csv(pao1_operon_filename, index_col=0, header=0)
# -

pao1_operons.head()

pao1_operons["Genes"] = pao1_operons["Genes"].str.split(";").apply(set)
pao1_operons.head()

# Remove operons with a single gene
pao1_operons = pao1_operons[pao1_operons["Genes"].apply(len) > 1]

gene_to_operons_df = pd.DataFrame(
    index=pao1_gene_module_labels.index, columns=list(pao1_operons.index)
)

# %%time
for gene in gene_to_operons_df.index:
    gene_to_operons_df.loc[gene] = [
        gene in pao1_operons.loc[operon, "Genes"] for operon in pao1_operons.index
    ]

# Add operons to pathway annotations
pao1_gene_annot = gene_to_pathways_df.merge(
    gene_to_operons_df, left_index=True, right_index=True, how="outer"
)

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

# ## Add regulon
#
# For each regulon, what genes are contained in it

# +
pao1_regulon_filename = "https://raw.githubusercontent.com/greenelab/core-accessory-interactome/6635c0e357c0172c2cebd0368648030e0ee4beaf/data/metadata/regulons_format.csv"

pao1_regulons = pd.read_csv(pao1_regulon_filename, index_col=0, header=0)
# -

pao1_regulons["Genes"] = pao1_regulons["Genes"].str.split(";").apply(set)

gene_to_regulons_df = pd.DataFrame(
    index=pao1_gene_module_labels.index, columns=list(pao1_regulons.index)
)

# %%time
for gene in gene_to_regulons_df.index:
    gene_to_regulons_df.loc[gene] = [
        gene in pao1_regulons.loc[regulon, "Genes"] for regulon in pao1_regulons.index
    ]

# Add regulons to other annotations
pao1_gene_annot = pao1_gene_annot.merge(
    gene_to_regulons_df, left_index=True, right_index=True, how="outer"
)

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

# ## Map pathway, operon, regulon to PA14
#
# The annotations we have are only for PAO1 genes, so we will map PAO1 core genes to PA14 core genes to add annotations to PA14

pao1_annotation_filename = paths.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(pao1_annotation_filename, "pao1")

gene_mapping_pao1 = gene_mapping_pao1["PA14_ID"].to_frame()

# Map PA14 gene ids
pao1_pa14_gene_annot = pao1_gene_annot.merge(
    gene_mapping_pao1, left_index=True, right_index=True
)

pao1_pa14_gene_annot.head()

# Reset index to PA14 gene ids
pa14_gene_annot = pao1_pa14_gene_annot.set_index("PA14_ID")
print(pa14_gene_annot.shape)
pa14_gene_annot.head()

# Merge annotations with module labels
pao1_gene_summary = pao1_gene_module_labels.merge(
    pao1_gene_annot, left_index=True, right_index=True, how="left"
)
pa14_gene_summary = pa14_gene_module_labels.merge(
    pa14_gene_annot, left_index=True, right_index=True, how="left"
)

print(pao1_gene_summary.shape)
pao1_gene_summary.head()

print(pa14_gene_summary.shape)
pa14_gene_summary.head()

# Drop duplicates
pa14_gene_summary = pa14_gene_summary[~pa14_gene_summary.index.duplicated(keep=False)]

# Save
pao1_gene_summary.to_csv(
    os.path.join(paths.LOCAL_DATA_DIR, f"pao1_gene_module_annotated_{method}.tsv"),
    sep="\t",
)
pa14_gene_summary.to_csv(
    os.path.join(paths.LOCAL_DATA_DIR, f"pa14_gene_module_annotated_{method}.tsv"),
    sep="\t",
)
