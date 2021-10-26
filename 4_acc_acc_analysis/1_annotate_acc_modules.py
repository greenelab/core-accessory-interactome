# -*- coding: utf-8 -*-
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
# 2. Baseline expression and expression in some context of interest
# 3. How clustered the module is on the genome
# 4. KEGG pathways that genes are found in
# 5. GO pathways genes are found in
# 6. Regulon/operon genes are found in
#
# All this information will help _P. aeruginosa_ experiments filter and determine which module might be interesting to explore.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import scipy.stats
import statsmodels.stats.multitest
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import paths, utils, modules, annotations

random.seed(1)
# -

# Clustering method used to obtain gene-module assignments
method = "affinity"
processed = "raw"

# +
# Import gene memberships
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method}_acc_{processed}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{method}_acc_{processed}.tsv"
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

# Import metadata of samples
metadata_filename = paths.SAMPLE_METADATA

# Get df with gene ids as indices and gene names as a column
# Having the data in a df instead of a series will just allow me to do my merges that are in the notebook
pao1_gene_annot = pao1_gene_annot["Name"].to_frame("gene name")
pa14_gene_annot = pa14_gene_annot["Name"].to_frame("gene name")

print(pao1_gene_annot.shape)
pao1_gene_annot.tail()

# +
# Use correlation matrix to get length of the genome
# -


# ## Add gene names
#
# **TO DO: Remove these PAO1 gene ids from the correlation and redo correlation and acc-acc. Then add how=left**

# Add gene names
pao1_gene_module_labels = pao1_membership.merge(
    pao1_gene_annot, left_index=True, right_index=True
)
pa14_gene_module_labels = pa14_membership.merge(
    pa14_gene_annot, left_index=True, right_index=True
)

# Note: Many gene ids don't have an associated gene name and so are NaNs
print(pao1_gene_module_labels.shape)
pao1_gene_module_labels.head()

# Note: Many gene ids don't have an associated gene name and so are NaNs
print(pa14_gene_module_labels.shape)
pa14_gene_module_labels.head()

# ## Add expression information
#
# 1. What is the baseline level of expression for each gene in the module?
# 2. What is the expression level of genes in a clinical context (i.e. clinical samples)?

# Read in expression data
# Data is of the form SRA sample id x gene id
pao1_compendium = pd.read_csv(paths.PAO1_COMPENDIUM, sep="\t", index_col=0)
pa14_compendium = pd.read_csv(paths.PA14_COMPENDIUM, sep="\t", index_col=0)

print(pao1_compendium.shape)
pao1_compendium.head()

print(pa14_compendium.shape)
pa14_compendium.head()

# Calculate median expression across all samples
pao1_median_all = pao1_compendium.median().to_frame("median expression")
pa14_median_all = pa14_compendium.median().to_frame("median expression")

pao1_median_all.head()

# +
# TO DO: Have Deb or Georgia select a study
# The following code blocks allow me to Select subset of samples and calculate the median
# expression across that subset of samples.
# An interesting selection would be what the clinical expression is, however
# it looks like we removed many of the clinical isolates from this compendium with our strain binning
# For now I will leave these blocks commented out
# selected_sample_ids = utils.get_sample_ids(
#   metadata_filename, experiment_colname="SRA_study", sample_colname="Experiment", experiment_id="SRP063289")

# +
# Subset compendium
# subset_pao1_compendium = pao1_compendium.loc[selected_sample_ids]
# subset_pa14_compendium = pa14_compendium.loc[selected_sample_ids]

# +
# print(subset_pao1_compendium.shape)
# print(subset_pa14_compendium.shape)

# +
# pao1_median_subset = subset_pao1_compendium.median().to_frame("median subset expression")
# pa14_median_subset = subset_pa14_compendium.median().to_frame("median subset expression")
# -

# Add median expression to gene ids
pao1_gene_annot = pao1_gene_module_labels.merge(
    pao1_median_all, left_index=True, right_index=True, how="left"
)
pa14_gene_annot = pa14_gene_module_labels.merge(
    pa14_median_all, left_index=True, right_index=True, how="left"
)

# Add median subset expression to gene ids
"""pao1_gene_annot = pao1_gene_annot.merge(
    pao1_median_subset, left_index=True, right_index=True, how="left"
)
pa14_gene_annot = pa14_gene_annot.merge(
    pa14_median_subset, left_index=True, right_index=True, how="left"
)"""

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

print(pa14_gene_annot.shape)
pa14_gene_annot.head()


# ## Genome location information
#
# How far are genes from other genes in the same module?

# +
# Sort gene ids and get last gene id to use as length of the genome
# This gene id should match the number of gene ids
sorted_pao1_compendium = pao1_compendium.T.sort_index()
pao1_last_gene_id = sorted_pao1_compendium.index[-1]

sorted_pa14_compendium = pa14_compendium.T.sort_index()
pa14_last_gene_id = sorted_pa14_compendium.index[-1]
# -

# Remove "PA" at the beginning of the identifier and convert into a float
pao1_genome_len = float(pao1_last_gene_id.split("PA")[-1])
pa14_genome_len = float(pa14_last_gene_id.split("PA14_")[-1])

print(pao1_genome_len, pa14_genome_len)

pao1_module_dist = modules.get_intra_module_dist(pao1_gene_annot, "PA", pao1_genome_len)
pa14_module_dist = modules.get_intra_module_dist(
    pa14_gene_annot, "PA14_", pa14_genome_len
)

pao1_module_dist.head(10)

pa14_module_dist.head(10)

# Add module distance to gene names
pao1_gene_annot = pao1_gene_annot.merge(
    pao1_module_dist, left_index=True, right_index=True, how="left"
)
pa14_gene_annot = pa14_gene_annot.merge(
    pa14_module_dist, left_index=True, right_index=True, how="left"
)

pao1_gene_annot.head()

# ## Add KEGG pathway enrichment analysis
#
# For each pathway, find significant association of pathways in accessory-accessory modules. This information is only available for PAO1.
#
# The [Fisher's exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) determines whether there is a significant association between two categorical variables in a contingency table (i.e two classifications of the data). Here we used use the Fisher’s exact test to determine if there is an association between the two classifications: in kegg pathway or not and in accessory-accessory module or not. In other words, we want to determine if there is a statistically significant association between genes found in a given accessory-accessory moudle and the genes involved in a given KEGG pathway. To do this we compare the ratio of genes found in the kegg pathway that are in the accessory-accessory module to the ratio of kegg pathway genes that are not found in the accessory-accessory module.
#
# Since the numbers are large, we also applied the $\chi^2$ test as an alternative to the Fisher's exact test.

# +
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

# pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", index_col=0, header=None)
# -

# pao1_pathways[2] = pao1_pathways[2].str.split(";").apply(set)
# pao1_pathways.index = pao1_pathways.index.str.split(" - ").str[0]
pao1_pathways = annotations.load_format_KEGG(pao1_pathway_filename)
pao1_pathways.head()

pao1_gene_annot.head()


# Given an accessory-accessory module, look for the array module with the most overlap/significant p-value
def KEGG_enrichment(acc_membership_df, kegg_df):
    all_genes = set(acc_membership_df.index)

    rows = []
    best_rows = []
    # For each accessory-accessory module
    for module_name, module_df_group in acc_membership_df.groupby("module id"):
        num_module_genes = module_df_group.shape[0]
        module_genes = set(module_df_group.index)
        not_module_genes = all_genes.difference(module_genes)

        # Find the KEGG pathway with the best overlap
        for kegg_name in kegg_df.index:
            num_kegg_genes = kegg_df.loc[kegg_name, 1]
            kegg_genes = set(kegg_df.loc[kegg_name, 2])
            not_kegg_genes = all_genes.difference(kegg_genes)

            # Make contingency table
            # -----------------|accessory module |not accessory module
            # kegg pathway     | # genes         | # genes
            # not kegg pathway | # genes         | # genes
            module_kegg_genes = module_genes.intersection(kegg_genes)
            not_module_kegg_genes = not_module_genes.intersection(kegg_genes)
            module_not_kegg_genes = module_genes.intersection(not_kegg_genes)
            not_module_not_kegg_genes = not_module_genes.intersection(not_kegg_genes)

            observed_contingency_table = np.array(
                [
                    [len(module_kegg_genes), len(not_module_kegg_genes)],
                    [len(module_not_kegg_genes), len(not_module_not_kegg_genes)],
                ]
            )

            # Fisher's exact test
            oddsr, pval = scipy.stats.fisher_exact(
                observed_contingency_table, alternative="greater"
            )
            # chi2 test will not accept 0 counts for the contingency table
            # chi2, pval, dof, expected_counts = scipy.stats.chi2_contingency(
            #    observed_contingency_table
            # )
            # print(oddsr, pval)

            rows.append(
                {
                    "module id": module_name,
                    "enriched KEGG pathway": kegg_name,
                    "p-value": pval,
                    "num shared genes": len(module_kegg_genes),
                    "size module": num_module_genes,
                    "size KEGG pathway": num_kegg_genes,
                }
            )

    enrichment_df = pd.DataFrame(rows)

    # Get corrected pvalues
    (
        reject_,
        pvals_corrected_,
        alphacSidak,
        alphacBonf,
    ) = statsmodels.stats.multitest.multipletests(
        enrichment_df["p-value"].values,
        alpha=0.05,
        method="fdr_bh",
        is_sorted=False,
    )

    enrichment_df["corrected p-value"] = pvals_corrected_

    # Select best module mapping
    for grp, grp_df in enrichment_df.groupby("module id"):
        # Find if any pathways is significant
        any_significant = (grp_df["corrected p-value"] < 0.05).any()
        if any_significant:
            best_kegg = grp_df[grp_df["corrected p-value"] < 0.05][
                "enriched KEGG pathway"
            ]
            best_pval = grp_df[grp_df["corrected p-value"] < 0.05]["p-value"].values[0]
            best_shared = grp_df[grp_df["corrected p-value"] < 0.05][
                "num shared genes"
            ].values[0]
            best_module_size = grp_df[grp_df["corrected p-value"] < 0.05][
                "size module"
            ].values[0]
            best_kegg_size = grp_df[grp_df["corrected p-value"] < 0.05][
                "size KEGG pathway"
            ].values[0]
            best_corrected_pval = grp_df[grp_df["corrected p-value"] < 0.05][
                "corrected p-value"
            ].values[0]
            best_rows.append(
                {
                    "module id": grp,
                    "enriched KEGG pathway": best_kegg,
                    "p-value": best_pval,
                    "num shared genes": best_shared,
                    "size module": best_module_size,
                    "size KEGG pathway": best_kegg_size,
                    "corrected p-value": best_corrected_pval,
                }
            )
        else:
            best_rows.append(
                {
                    "module id": grp,
                    "enriched KEGG pathway": "NA",
                    "p-value": "NA",
                    "num shared genes": "NA",
                    "size module": "NA",
                    "size KEGG pathway": "NA",
                    "corrected p-value": "NA",
                }
            )
    best_enrichment_df = pd.DataFrame(best_rows).set_index("module id")

    return best_enrichment_df


pao1_enrichment_df = KEGG_enrichment(pao1_membership, pao1_pathways)

pao1_enrichment_df.head(20)

# Add pathway enrichment information
pao1_gene_annot = pao1_gene_annot.merge(
    pao1_enrichment_df, left_on="module id", right_index=True, how="left"
)

pao1_gene_annot.head()

# ## Import and format operon

pao1_operon_filename = paths.PAO1_OPERON
pa14_operon_filename = paths.PA14_OPERON

pao1_operon = annotations.load_format_operons(pao1_operon_filename)
pa14_operon = annotations.load_format_operons(pa14_operon_filename)

pao1_operon.head()

# Add operons to pathway annotations for PAO1
pao1_gene_annot = pao1_gene_annot.merge(
    pao1_operon, left_index=True, right_index=True, how="left"
)

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

# For PA14 we only have operon annotations
pa14_gene_annot = pa14_gene_annot.merge(
    pa14_operon, left_index=True, right_index=True, how="left"
)

# ## Add regulon
#
# For each regulon, what genes are contained in it. This information is only available for PAO1

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
    gene_to_regulons_df, left_index=True, right_index=True, how="left"
)

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

print(pa14_gene_annot.shape)
pa14_gene_annot.head()

# ## Plot trends

# +
# Get pairwise distance per module
# Otherwise plotting will be scaled by the number of genes in a module
pao1_mean_dist = []
for grp_name, grp_df in pao1_gene_annot.groupby("module id"):
    pao1_mean_dist.append(grp_df["mean pairwise dist"][0])

pa14_mean_dist = []
for grp_name, grp_df in pa14_gene_annot.groupby("module id"):
    pa14_mean_dist.append(grp_df["mean pairwise dist"][0])
# -

sns.displot(pao1_mean_dist)
plt.title("Mean pairwise distance of PAO1 accessory modules")

sns.displot(pa14_mean_dist)
plt.title("Mean pairwise distance of PA14 accessory modules")

# Find which modules/genes have a large range
pao1_farapart_genes = pao1_gene_annot[pao1_gene_annot["mean pairwise dist"] >= 2000]

pa14_farapart_genes = pa14_gene_annot[pa14_gene_annot["mean pairwise dist"] >= 20000]

# Save
pao1_gene_annot.to_csv(f"pao1_acc_gene_module_annotated_{method}.tsv", sep="\t")
pa14_gene_annot.to_csv(f"pa14_acc_gene_module_annotated_{method}.tsv", sep="\t")
pao1_farapart_genes.to_csv(f"pao1_farapart_acc_modules_{method}.tsv", sep="\t")
pa14_farapart_genes.to_csv(f"pa14_farapart_acc_modules_{method}.tsv", sep="\t")

# These annotations will be used to help _P. aeruginosa_ experts, like our collaborators, to determine what accessory-accessory modules to focus on.
#
#
# Note: Since genes can be in multiple KEGG pathways and regulons, each pathway and regulon are separate columns. Whereas operons are a single column since genes can belong to only a single operon.
#
# Maybe focus on the ones with a high range since those might tell us about how co-regulation evolves. These far apart come to be regulated by the same one vs all genes that evolved together?
