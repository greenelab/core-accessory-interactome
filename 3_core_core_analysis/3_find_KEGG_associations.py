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

# # Find KEGG associations
#
# This notebook will create a table that has the KEGG pathways that are associated with all genes, but we are particularly interested in those that are associated with the most and least stable genes.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
from scripts import paths, utils, modules, annotations

random.seed(1)
# -

# Output files
pao1_out_filename = "pao1_core_similarity_associations_spell.tsv"
pa14_out_filename = "pa14_core_similarity_associations_spell.tsv"

# +
# Load transcriptional similarity df
pao1_similarity_scores_filename = "pao1_core_similarity_expression_stats_spell.tsv"
pa14_similarity_scores_filename = "pa14_core_similarity_expression_stats_spell.tsv"

pao1_similarity_scores = pd.read_csv(
    pao1_similarity_scores_filename, sep="\t", header=0, index_col=0
)
pa14_similarity_scores = pd.read_csv(
    pa14_similarity_scores_filename, sep="\t", header=0, index_col=0
)
# -

print(pao1_similarity_scores.shape)
pao1_similarity_scores.head()

print(pa14_similarity_scores.shape)
pa14_similarity_scores.head()

# ## Load KEGG annotations

pao1_pathway_filename = "pao1_kegg_annot.tsv"
pa14_pathway_filename = "pa14_kegg_annot.tsv"

pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", header=0, index_col=0)
pa14_pathways = pd.read_csv(pa14_pathway_filename, sep="\t", header=0, index_col=0)

print(pao1_pathways.shape)
pao1_pathways.head()

print(pa14_pathways.shape)
pa14_pathways.head()


# ## Get pathway associations for all genes

def get_associated_pathways(genes_, pathway_df):
    rows = []
    for gene_id in genes_:
        pathway_bool = [
            gene_id in pathway_df.loc[pathway, "gene_ids"]
            for pathway in pathway_df.index
        ]
        found_pathways = list(pathway_df[pathway_bool].index)
        rows.append({"gene id": gene_id, "pathways present": found_pathways})
    return pd.DataFrame(rows).set_index("gene id")


# Get KEGG associations for all genes in PAO1
all_pao1_gene_ids = list(pao1_similarity_scores.index)
pao1_associations = get_associated_pathways(all_pao1_gene_ids, pao1_pathways)

print(pao1_associations.shape)
pao1_associations.head()

# Get KEGG associations for all genes in PA14
all_pa14_gene_ids = list(pa14_similarity_scores.index)
pa14_associations = get_associated_pathways(all_pa14_gene_ids, pa14_pathways)

print(pa14_associations.shape)
pa14_associations.head()

# Merge KEGG associations with transcriptional similarity information
pao1_associations = pao1_similarity_scores.merge(
    pao1_associations, left_index=True, right_index=True, how="left"
)
pa14_associations = pa14_similarity_scores.merge(
    pa14_associations, left_index=True, right_index=True, how="left"
)

print(pao1_associations.shape)
pao1_associations.head()

print(pa14_associations.shape)
pa14_associations.head()

# Save
pao1_associations.to_csv(pao1_out_filename, sep="\t")
pa14_associations.to_csv(pa14_out_filename, sep="\t")
