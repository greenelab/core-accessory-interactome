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

# # Common stability analysis
#
# This notebook examines the relationship between commonly DE genes and the stability of those genes

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import utils, paths

# ### Get common DE statistics

# +
# Load summary statistics generated by SOPHIE using the PAO1 and PA
pao1_common_DEGs_filename = "find_common_DEGs/generic_gene_summary_SRP117105.tsv"
pa14_common_DEGs_filename = "find_common_DEGs/generic_gene_summary_SRP074292.tsv"


pao1_SOPHIE_stats = pd.read_csv(
    pao1_common_DEGs_filename, sep="\t", index_col=0, header=0
)
pa14_SOPHIE_stats = pd.read_csv(
    pa14_common_DEGs_filename, sep="\t", index_col=0, header=0
)
# -

print(pao1_SOPHIE_stats.shape)
pao1_SOPHIE_stats.head()

# ### Get stability statistics

# +
# Load transcriptional similarity df
# These are the subset of genes that we will consider
pao1_similarity_scores_filename = "../3_core_core_analysis/pao1_similarity_scores.tsv"
pa14_similarity_scores_filename = "../3_core_core_analysis/pa14_similarity_scores.tsv"

pao1_similarity_scores = pd.read_csv(
    pao1_similarity_scores_filename, sep="\t", header=0, index_col=0
)
pa14_similarity_scores = pd.read_csv(
    pa14_similarity_scores_filename, sep="\t", header=0, index_col=0
)
# -

print(pao1_similarity_scores.shape)
pao1_similarity_scores.head()

# Merge transcriptional stability statistics and common statistics
pao1_all_stats = pao1_similarity_scores.merge(
    pao1_SOPHIE_stats,
    left_index=True,
    right_index=True,
)
pa14_all_stats = pa14_similarity_scores.merge(
    pa14_SOPHIE_stats,
    left_index=True,
    right_index=True,
)

print(pao1_all_stats.shape)
pao1_all_stats.head()

sns.jointplot(
    data=pao1_all_stats,
    x="Transcriptional similarity across strains",
    y="Percentile (simulated)",
    kind="hex",
    # alpha=0.2
)
plt.suptitle("PAO1 core gene stability vs commonality", y=1.05)

sns.jointplot(
    data=pao1_all_stats,
    x="Transcriptional similarity across strains",
    y="Z score",
    kind="hex",
)
plt.suptitle("PAO1 core gene stability vs commonality (z-score)", y=1.05)

sns.jointplot(
    data=pa14_all_stats,
    x="Transcriptional similarity across strains",
    y="Percentile (simulated)",
    kind="hex",
)
plt.suptitle("PA14 core gene stability vs commonality", y=1.05)
