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

# # Decide threshold
#
# The goal of this notebook is to define a threshold to determine if a sample if PAO1 or not (likewise, if a sample is PA14 or not). We used known labels from SRA to do this. Specifically, we examined the distribution of PAO1 samples (grey) vs non-PAO1 samples (blue). We define the threshold to be one that separated between the two distributions. We use this threshold in [2_create_compendia.ipynb](2_create_compendia.ipynb) to partition gene expression data into PAO1 and PA14 compendia.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import seaborn as sns
from textwrap import fill
import matplotlib.pyplot as plt
from scripts import paths, utils

# +
# Raw (normalized counts) expression data files
pao1_expression_filename = paths.PAO1_GE
pa14_expression_filename = paths.PA14_GE

# File containing table to map sample id to strain name
sample_to_strain_filename = paths.SAMPLE_TO_STRAIN
# -

# Load expression data
pao1_expression = pd.read_csv(pao1_expression_filename, sep="\t", index_col=0, header=0)
pa14_expression = pd.read_csv(pa14_expression_filename, sep="\t", index_col=0, header=0)

# Load metadata
# Set index to experiment id, which is what we will use to map to expression data
sample_to_strain_table_full = pd.read_csv(sample_to_strain_filename, index_col=2)

# +
threshold = 25

pao1_dist_filename = (
    f"dist_median_acc_expression_pao1_compendium_{threshold}threshold.svg"
)
pa14_dist_filename = (
    f"dist_median_acc_expression_pa14_compendium_{threshold}threshold.svg"
)
# -

# ## Format SRA annotations

# +
# Since experiments have multiple runs there are duplicated experiment ids in the index
# We will need to remove these so that the count calculations are accurate
sample_to_strain_table_full_processed = sample_to_strain_table_full[
    ~sample_to_strain_table_full.index.duplicated(keep="first")
]

assert (
    len(sample_to_strain_table_full.index.unique())
    == sample_to_strain_table_full_processed.shape[0]
)

# +
# Aggregate boolean labels into a single strain label
aggregated_label = []
for exp_id in list(sample_to_strain_table_full_processed.index):
    if sample_to_strain_table_full_processed.loc[exp_id, "PAO1"].all() == True:
        aggregated_label.append("PAO1")
    elif sample_to_strain_table_full_processed.loc[exp_id, "PA14"].all() == True:
        aggregated_label.append("PA14")
    elif sample_to_strain_table_full_processed.loc[exp_id, "PAK"].all() == True:
        aggregated_label.append("PAK")
    elif (
        sample_to_strain_table_full_processed.loc[exp_id, "ClinicalIsolate"].all()
        == True
    ):
        aggregated_label.append("Clinical Isolate")
    else:
        aggregated_label.append("NA")

sample_to_strain_table_full_processed["Strain type"] = aggregated_label

sample_to_strain_table = sample_to_strain_table_full_processed["Strain type"].to_frame()

sample_to_strain_table.head()
# -

# ## Label samples with SRA annotations

# Label samples with SRA annotations
# pao1_expression_label = pao1_expression_binned.join(
#    sample_to_strain_table, how='left')
pao1_expression_label = pao1_expression.merge(
    sample_to_strain_table, left_index=True, right_index=True
)
pa14_expression_label = pa14_expression.merge(
    sample_to_strain_table, left_index=True, right_index=True
)
print(pao1_expression_label.shape)
pao1_expression_label.head()

print(pa14_expression_label.shape)
pa14_expression_label.head()

# ## Get accessory gene expression

# +
pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT

core_acc_dict = utils.get_my_core_acc_genes(
    pao1_annot_filename, pa14_annot_filename, pao1_expression, pa14_expression
)
# -

pao1_acc = core_acc_dict["acc_pao1"]
pa14_acc = core_acc_dict["acc_pa14"]

# +
# Create accessory df
# accessory gene ids | median accessory expression | strain label

# PAO1
pao1_acc_expression = pao1_expression_label[pao1_acc]
pao1_acc_expression["median_acc_expression"] = pao1_acc_expression.median(axis=1)

# PA14
pa14_acc_expression = pa14_expression_label[pa14_acc]
pa14_acc_expression["median_acc_expression"] = pa14_acc_expression.median(axis=1)
# -

# Add back labels
pao1_acc_expression["Strain type"] = pao1_expression_label["Strain type"]
pa14_acc_expression["Strain type"] = pa14_expression_label["Strain type"]

print(pao1_acc_expression.shape)
pao1_acc_expression.head()

print(pa14_acc_expression.shape)
pa14_acc_expression.head()

# +
# Merge PAO1 and PA14 accessory dataframes
pao1_pa14_acc_expression = pao1_acc_expression.merge(
    pa14_acc_expression,
    left_index=True,
    right_index=True,
    suffixes=["_pao1", "_pa14"],
)

print(pao1_pa14_acc_expression.shape)
pao1_pa14_acc_expression.head()
# -

# ## Plot distribution

# +
# Get PAO1 samples that are labeled PAO1 and non-PAO1
pao1_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pao1"] == "PAO1",
    "median_acc_expression_pao1",
]

non_pao1_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pao1"] != "PAO1",
    "median_acc_expression_pao1",
]
# -

print(pao1_sra.shape)
pao1_sra.head()

print(non_pao1_sra.shape)
non_pao1_sra.head()

646 + 1687

# +
f = sns.distplot(
    non_pao1_sra,
    color="#795C34",
    label="non-PAO1",
    kde=False,
    hist_kws={
        "alpha": 0.6,
    },
)
f = sns.distplot(
    pao1_sra, color="#C6A9B5", label="PAO1", kde=False, hist_kws={"alpha": 0.7}
)
plt.axvline(threshold, color="black", linestyle="--")
f.set_yscale("log")
f.set_ylabel("Count", fontsize=18)
f.set_xlabel("PAO1 expression", fontsize=18)
f.tick_params(labelsize=16)
plt.legend(fontsize=16)

f.figure.savefig(pao1_dist_filename, bbox_inches="tight", format="svg", dpi=300)

# +
# Get PA14 samples that are labeled PA14 and non-PA14
pa14_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pa14"] == "PA14",
    "median_acc_expression_pa14",
]

non_pa14_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pa14"] != "PA14",
    "median_acc_expression_pa14",
]
# -

print(pa14_sra.shape)
pa14_sra.head()

print(non_pa14_sra.shape)
non_pa14_sra.head()

441 + 1892

# +
g = sns.distplot(
    non_pa14_sra,
    color="#795C34",
    label="non-PA14",
    kde=False,
    hist_kws={"alpha": 0.6},
)
g = sns.distplot(
    pa14_sra,
    color="#895881",
    label="PA14",
    kde=False,
    hist_kws={"alpha": 0.7},
)

plt.axvline(threshold, color="black", linestyle="--")
g.set_yscale("log")
g.set_ylabel("Count", fontsize=18)
g.set_xlabel("PA14 expression", fontsize=18)
g.tick_params(labelsize=16)
plt.legend(fontsize=16)
g.figure.savefig(pa14_dist_filename, bbox_inches="tight", format="svg", dpi=300)
# -

# **Takeaway:**
# Looks like using a threshold of 25 normalized counts separates between SRA-annotated PAO1 samples vs non-PAO1 samples. Similarly for PA14. This is the threshold we'll use to bin samples into PAO1 vs PA14 compendia.

# +
# Save df with median accessory gene expression for user resource
pao1_pa14_acc_expression_out = pao1_pa14_acc_expression.drop(
    columns=["Strain type_pao1", "Strain type_pa14"]
)

pao1_pa14_acc_expression_out.to_csv("median_acc_expression.tsv", sep="\t")
