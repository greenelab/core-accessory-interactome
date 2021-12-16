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

# # Examine samples with high PAO1 and PA14 accessory expression
#
# Based on the scatterplot in the [previous notebook](cluster_by_accessory_gene.ipynb), there are some samples with high PAO1 accessory and high PA14 accessory expression. We want to know what these samples are

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
from scripts import paths, utils

# Expression data pre-binning
pao1_expression_prebin_filename = paths.PAO1_PREBIN_COMPENDIUM
pa14_expression_prebin_filename = paths.PA14_PREBIN_COMPENDIUM

# Load data
pao1_expression_prebin = pd.read_csv(
    pao1_expression_prebin_filename, sep="\t", index_col=0, header=0
)
pa14_expression_prebin = pd.read_csv(
    pa14_expression_prebin_filename, sep="\t", index_col=0, header=0
)

# ## Get core/accessory genes

# +
pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT

core_acc_dict = utils.get_my_core_acc_genes(
    pao1_annot_filename,
    pa14_annot_filename,
    pao1_expression_prebin,
    pa14_expression_prebin,
)
# -

pao1_core = core_acc_dict["core_pao1"]
pa14_core = core_acc_dict["core_pa14"]
pao1_acc = core_acc_dict["acc_pao1"]
pa14_acc = core_acc_dict["acc_pa14"]

# ## Find accessory gene expression

# +
# Create accessory df for PAO1 compendium
# accessory gene ids | median accessory expression | strain label

# PAO1-only genes in PAO1 compendium
pao1_acc_pao1_compendium = pao1_expression_prebin[pao1_acc]
pao1_acc_pao1_compendium["median acc expression"] = pao1_acc_pao1_compendium.median(
    axis=1
)

# +
# Create accessory df for PA14 compendium
# accessory gene ids | median accessory expression | strain label

# PA14-only genes in PA14 compendium
pa14_acc_pa14_compendium = pa14_expression_prebin[pa14_acc]
pa14_acc_pa14_compendium["median acc expression"] = pa14_acc_pa14_compendium.median(
    axis=1
)
# -

high_pao1_acc_samples = pao1_acc_pao1_compendium[
    pao1_acc_pao1_compendium["median acc expression"] > 40
].index
high_pa14_acc_samples = pa14_acc_pa14_compendium[
    pa14_acc_pa14_compendium["median acc expression"] > 40
].index

shared_samples = set(high_pao1_acc_samples).intersection(high_pa14_acc_samples)

print(len(shared_samples))
shared_samples

# Who are these samples?
# * [SRX6437734](https://www.ncbi.nlm.nih.gov/sra/?term=SRX6437734) and SRX6437785: From an expteriment comparing transcriptional patterns in biofilm vs planktonic growth. This sample is a biofilm sample.
# * [SRX6911057](https://www.ncbi.nlm.nih.gov/sra/?term=SRX6911057) a clinical isolate.
#
# There are not too many that had a high expression in both but 2 are from the same experiment studying biofilm growth and the 1 is a clinical isolate.
