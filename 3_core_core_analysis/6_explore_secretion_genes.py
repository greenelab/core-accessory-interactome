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

# # Explore secretion system genes
#
# [KEGG enrichment analysis](5_KEGG_enrichment_of_stable_genes.ipynb) found that genes associated with ribosome, Lipopolysaccharide (outer membrane) biosynthesis, citrate cycle are significantly conserved across strains.
# Indeed functions that are essential seem to be significantly conserved across strains as expected.
# However, there are also pathways like the secretion systems, which allow for inter-strain warfare, that we’d expect to vary across strains but were found to be conserved (T3SS significant but not T6SS).
#
# This notebook examines the stability score of the genes in the secretion systems to determine if there is a subset of the secretion genes, related to the machinery that is conserved while others, like the secretory proteins, are more variable.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import annotations

random.seed(1)
# -

# ## Load data and metadata

# Input similarity scores and annotations filenames
# Since the results are similar we only need to look at the scores for one strain type
pao1_similarity_filename = "pao1_core_similarity_associations_final_spell.tsv"

# Import df
pao1_similarity = pd.read_csv(pao1_similarity_filename, sep="\t", index_col=0, header=0)

pao1_similarity.head()

# Load KEGG pathway data
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

pao1_pathways = annotations.load_format_KEGG(pao1_pathway_filename)
print(pao1_pathways.shape)
pao1_pathways.head()

# ## Get genes related to secretion pathways

pao1_pathways.loc[
    [
        "KEGG-Module-M00334: Type VI secretion system",
        "KEGG-Module-M00332: Type III secretion system",
        "KEGG-Module-M00335: Sec (secretion) system",
    ]
]

# Get genes related to pathways
T6SS_genes = list(pao1_pathways.loc["KEGG-Module-M00334: Type VI secretion system", 2])
T3SS_genes = list(pao1_pathways.loc["KEGG-Module-M00332: Type III secretion system", 2])
secretion_genes = list(
    pao1_pathways.loc["KEGG-Module-M00335: Sec (secretion) system", 2]
)

# Pull out genes related to T3SS
T6SS_similarity = pao1_similarity.reindex(T6SS_genes)
T3SS_similarity = pao1_similarity.reindex(T3SS_genes)
sec_similarity = pao1_similarity.reindex(secretion_genes)

T6SS_similarity.sort_values(by="Transcriptional similarity across strains")

T3SS_similarity.sort_values(by="Transcriptional similarity across strains")

# +
# sec_similarity.sort_values(by="Transcriptional similarity across strains")
# -

# Save T3SS and T6SS df for easier lookup
T3SS_similarity.to_csv("T3SS_core_similarity_associations_final_spell.tsv", sep="\t")
T6SS_similarity.to_csv("T6SS_core_similarity_associations_final_spell.tsv", sep="\t")

# ## Plot

# +
plt.figure(figsize=(10, 8))
sns.violinplot(
    data=pao1_similarity,
    x="Transcriptional similarity across strains",
    palette="Blues",
    inner=None,
)
sns.swarmplot(
    data=T6SS_similarity,
    x="Transcriptional similarity across strains",
    color="k",
    label="T6SS genes",
    alpha=0.8,
)

sns.swarmplot(
    data=T3SS_similarity,
    x="Transcriptional similarity across strains",
    color="r",
    label="T3SS genes",
    alpha=0.8,
)

# sns.swarmplot(
#    data=sec_similarity,
#    x="Transcriptional similarity across strains",
#    color="yellow",
#    label="secretion system genes",
#    alpha=0.8,
# )

# Add text labels for least stable genes amongst the T3SS/T6SS
plt.text(
    x=T3SS_similarity.loc[
        T3SS_similarity["Name"] == "pscR", "Transcriptional similarity across strains"
    ],
    y=0.02,
    s="$pscR$",
)

plt.text(
    x=T6SS_similarity.loc[
        T6SS_similarity["Name"] == "vgrG6", "Transcriptional similarity across strains"
    ],
    y=-0.02,
    s="$vgrG6$",
)
plt.text(
    x=T6SS_similarity.loc[
        T6SS_similarity["Name"] == "vgrG3", "Transcriptional similarity across strains"
    ],
    y=-0.02,
    s="$vgrG3$",
)
plt.text(
    x=T6SS_similarity.loc[
        T6SS_similarity["Name"] == "vgrG4a", "Transcriptional similarity across strains"
    ],
    y=-0.02,
    s="$vgrG4a$",
)

plt.title("Stability of secretion system genes", fontsize=14)
plt.legend()
# -

# We hypothesized that most secretion machinery genes would be conserved but that secreted proteins (i.e. effector proteins) would be less conserved. In general, the effector proteins were not included in the KEGG annotations, which is probably why these secretion systems were found to be highly stable.
#
# T6SS genes are among the most stable with the vgrG6, vgrG3, vgrG4a among the least stable. T3SS among the most stable with pscR among the least stable
#
# Need to read more about these genes and if it makes sense that they are at the bottom.
