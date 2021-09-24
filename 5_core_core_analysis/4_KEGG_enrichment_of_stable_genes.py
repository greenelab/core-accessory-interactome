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
import scipy.stats
import statsmodels.stats.multitest
from scripts import paths, utils, annotations

# Load KEGG pathway data
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

pao1_pathways = annotations.load_format_KEGG(pao1_pathway_filename)
print(pao1_pathways.shape)
pao1_pathways.head()

# +
# Load transcriptional similarity df
# These are the subset of genes that we will consider
pao1_similarity_scores_filename = "../5_core_core_analysis/pao1_similarity_scores.tsv"

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
        num_kegg_genes = kegg_df.loc[kegg_name, 1]
        kegg_genes = set(kegg_df.loc[kegg_name, 2])
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

# Save
pao1_most_stable_enrichment.to_csv("pao1_most_stable_enrichment.tsv", sep="\t")
pao1_least_stable_enrichment.to_csv("pao1_least_stable_enrichment.tsv", sep="\t")