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

# # PAO1 vs PA14 core modules
#
# This notebook compares the core modules from PAO1 and PA14 and tries to determine if they are consistent (i.e. are the same genes grouped in PAO1 also grouped in PA14?)

# %load_ext autoreload
# %autoreload 2
import os
import pandas as pd
from core_acc_modules import paths, utils

# User params
cluster_method = "dbscan"

# Load module membership
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{cluster_method}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{cluster_method}.tsv"
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", index_col=0, header=0)
pa14_membership = pd.read_csv(pa14_membership_filename, sep="\t", index_col=0, header=0)

print(pao1_membership.shape)
pao1_membership.head()

print(pa14_membership.shape)
pa14_membership.head()

# ## Get gene id mapping
#
# Get mapping from PAO1 gene ids to PA14 gene ids and vice versa in order to map core modules from PAO1 to PA14 and vice versa

pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT

# Map PA14 gene ids to PAO1 gene ids
pao1_to_pa14_map = utils.get_pao1_pa14_gene_map(pao1_annot_filename, "PAO1")
pa14_to_pao1_map = utils.get_pao1_pa14_gene_map(pao1_annot_filename, "PA14")

pao1_to_pa14_map_dict = pao1_to_pa14_map["PA14_ID"].to_dict()
pa14_to_pao1_map_dict = pa14_to_pao1_map["PAO1_ID"].to_dict()

# Test
# pa14_to_pao1_map.loc["PA14_72210"]
pao1_to_pa14_map.loc["PA0052"]

# Assign genes in PAO1 modules to PA14 gene ids
pao1_membership["index mapped"] = (
    pao1_membership.reset_index()["index"].map(pao1_to_pa14_map_dict).values
)

# Assign genes in PA14 modules to PAO1 gene ids
pa14_membership["index mapped"] = (
    pa14_membership.reset_index()["index"].map(pa14_to_pao1_map_dict).values
)

print(pa14_membership.shape)
pa14_membership.head()

print(pao1_membership.shape)
pao1_membership.head()

# ## Find modules that best match

# +
# PAO1 --> PA14
num_pao1_genes = len(pao1_core)

# For each module in PAO1
for pao1_group_name, pao1_df_group in pao1_membership.groupby("module id"):
    print(pao1_group_name)
    print(pao1_df_group)

    # Find the PA14 module with the best overlap
    for pa14_group_name, pa14_df_group in pa14_membership.groupby("module id"):
        print(pa14_group_name)
        print(pa14_df_group)

        # TO DO
        # Update based this analysis
        # Save p-values and best group_name
        num_generic_Crow_genes = shared_ranking.query(f"{ref_rank_col}>=80.0").shape[0]
        num_generic_SOPHIE_genes = shared_ranking[
            shared_ranking["Percentile (simulated)"] >= percentile_threshold
        ].shape[0]
        num_concordant_generic_genes = shared_ranking[
            (shared_ranking[ref_rank_col] >= percentile_threshold)
            & (shared_ranking["Percentile (simulated)"] >= percentile_threshold)
        ].shape[0]

        print(num_Crow_genes)
        print(num_generic_Crow_genes)
        print(num_generic_SOPHIE_genes)
        print(num_concordant_generic_genes)

        p = ss.hypergeom.sf(
            num_concordant_generic_genes,
            num_pao1_genes,
            num_generic_Crow_genes,
            num_generic_SOPHIE_genes,
        )
        print(p)
        break
    break
    # Store p-values

# +
# TO DO
# Plot pvalue dist
