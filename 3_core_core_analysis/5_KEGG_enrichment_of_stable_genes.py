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

# # KEGG enrichment of stable genes
#
# This notebooks looks at the group of most and least stable genes and performs a KEGG enrichment analysis to determine if there are any KEGG pathways that are significantly over-represented in our most or least stable gene sets.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import json
import statsmodels.stats.multitest
from scripts import paths, utils, annotations

# Load KEGG pathway data
pao1_pathway_filename = "pao1_kegg_annot.tsv"

pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", header=0, index_col=0)

print(pao1_pathways.shape)
pao1_pathways.head()

# +
# Load transcriptional similarity df
# These are the subset of genes that we will consider
pao1_similarity_scores_filename = "pao1_similarity_scores_spell.tsv"

pao1_similarity_scores = pd.read_csv(
    pao1_similarity_scores_filename, sep="\t", header=0, index_col=0
)
# -

pao1_similarity_scores.head()

# Get most and least stable core genes based on label
pao1_most_stable_genes = list(
    pao1_similarity_scores[pao1_similarity_scores["label"] == "most stable"].index
)
pao1_least_stable_genes = list(
    pao1_similarity_scores[pao1_similarity_scores["label"] == "least stable"].index
)


# For each KEGG pathway, perform stat test, save p-values to get corrected p-values, report stats per pathway
def KEGG_enrichment_of_stable_genes(similarity_score_df, gene_list, kegg_df):
    """
    This function performs a KEGG enrichment using most or least stable genes,
    provided in `gene_list`
    """

    all_genes = set(similarity_score_df.index)
    module_genes = set(gene_list)
    not_module_genes = all_genes.difference(module_genes)

    rows = []
    # Find the KEGG pathway with significant over-representation
    for kegg_name in kegg_df.index:
        num_kegg_genes = kegg_df.loc[kegg_name, "num_genes"]
        kegg_genes = kegg_df.loc[kegg_name, "gene_ids"]
        kegg_genes = set(json.loads(kegg_genes.replace("'", '"')))

        not_kegg_genes = all_genes.difference(kegg_genes)

        # Make contingency table
        # ---------------------| most stable  | not most stable
        # in KEGG pathway      | # genes      | # genes
        # not in KEGG pathway  | # genes     | # genes
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
                "enriched KEGG pathway": kegg_name,
                "odds ratio": oddsr,
                "p-value": pval,
                "num shared genes": len(module_kegg_genes),
                "size gene set": len(module_genes),
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
    enrichment_df["-log10 adj p-value"] = -np.log10(enrichment_df["corrected p-value"])

    return enrichment_df


pao1_most_stable_enrichment = KEGG_enrichment_of_stable_genes(
    pao1_similarity_scores, pao1_most_stable_genes, pao1_pathways
)

pao1_least_stable_enrichment = KEGG_enrichment_of_stable_genes(
    pao1_similarity_scores, pao1_least_stable_genes, pao1_pathways
)

print(pao1_most_stable_enrichment.shape)
pao1_most_stable_enrichment.sort_values(by="corrected p-value").head()

print(pao1_least_stable_enrichment.shape)
pao1_least_stable_enrichment.sort_values(by="corrected p-value").head()

# ## Plot

# +
# Create a shared scale to use for both plots
# np.linspace(3.759057e-22, 1)

# +
# Plot top most stable
pao1_most_stable_enrichment_top = pao1_most_stable_enrichment.sort_values(
    by="corrected p-value"
).head(10)
norm = plt.Normalize(0, 22)
cmap = sns.cubehelix_palette(as_cmap=True)

plt.figure(figsize=(8, 6))
f = plt.scatter(
    data=pao1_most_stable_enrichment_top,
    x="odds ratio",
    y="enriched KEGG pathway",
    s=300,
    c="-log10 adj p-value",
    norm=norm,
    cmap=cmap,
)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cf = plt.colorbar(sm)
cf.ax.set_title("-log10 adj p-value", fontsize=12)

plt.title("Enrichment of most stable core genes", fontsize=16, y=1.05)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("Odds ratio", fontsize=16)
plt.ylabel("")

# +
# Plot top most stable
pao1_least_stable_enrichment_top = pao1_least_stable_enrichment.sort_values(
    by="corrected p-value"
).head(10)
plt.figure(figsize=(8, 6))
norm = plt.Normalize(0.0, 22)
cmap = sns.cubehelix_palette(as_cmap=True)

g = plt.scatter(
    data=pao1_least_stable_enrichment_top,
    x="odds ratio",
    y="enriched KEGG pathway",
    s=300,
    c="-log10 adj p-value",
    norm=norm,
    cmap=cmap,
)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cg = plt.colorbar(sm)
cg.ax.set_title("-log10 adj p-value", fontsize=12)

plt.title("Enrichment of least stable core genes", fontsize=16, y=1.05)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("Odds ratio", fontsize=16)
plt.ylabel("")

# TO DO
# Update colorbar text to include label for -log10 p-value
# -

# Save
pao1_most_stable_enrichment.to_csv("pao1_most_stable_enrichment_spell.tsv", sep="\t")
pao1_least_stable_enrichment.to_csv("pao1_least_stable_enrichment_spell.tsv", sep="\t")

# **Takeaway:**
#
# KEGG enrichment analysis found that stable genes were significantly associated with essential functions: ribosome, lipopolysaccharide biosynthesis, citrate cycle. While least stable core genes are enriched in pathways related to degredation of aromatic compounds
#
# KEGG enrichment: https://docs.google.com/spreadsheets/d/1lXZZXXjZSOuQ-cMOZ9I5llIJ2OBOQBur0Spso2WN4oY/edit#gid=0
