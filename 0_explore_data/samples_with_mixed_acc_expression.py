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
# There are some samples with high PAO1 accessory and high PA14 accessory expression. We want to know what these samples are

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
    pao1_acc_pao1_compendium["median acc expression"] > 50
].index
high_pa14_acc_samples = pa14_acc_pa14_compendium[
    pa14_acc_pa14_compendium["median acc expression"] > 50
].index

shared_samples = set(high_pao1_acc_samples).intersection(high_pa14_acc_samples)

print(len(shared_samples))
shared_samples

# Who are these samples?
# * 'ERX2813655' Negative control of RpoS synthesis, which is involved in QS: https://www.ncbi.nlm.nih.gov/sra/?term=ERX2813655
# * 2 samples ('SRX1096902', 'SRX1097013') from a study testing effect of antibiotic treatment. 'SRX1096902' is 0m after treatment, 'SRX1097013' is 30m after treatment: https://www.ncbi.nlm.nih.gov/sra/?term=SRX1096902; https://www.ncbi.nlm.nih.gov/sra/?term=SRX1097013
# * 'SRX1127442' sRNA content of P. aeruginosa OMVs compared to whole cells: https://www.ncbi.nlm.nih.gov/sra/?term=SRX1127442
# * 2 samples ('SRX1491235', 'SRX1491236') testing transcriptomic response of PAO1 under Elevated Temperature: https://www.ncbi.nlm.nih.gov/sra/?term=SRX1491235
# * Sample ('SRX1516058') from a study reports the RNA-Seq analysis of bacteriophage PAK_P3 infecting PAK strain of P. aeruginosa: https://www.ncbi.nlm.nih.gov/sra/?term=SRX1516058
# * Sample ('SRX1747042') characterzing the in vitro and in vivo transcriptome of PAO1
# * Sample ('SRX2662725', 'SRX2662726', 'SRX2662727', 'SRX2662728') characterizing the transcriptional patterns of PAO1 under different culture conditions: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA360055
# * 4 samples ('SRX5000019', 'SRX5000020', 'SRX5000025', 'SRX5000026') co-cultured with Aspergillus fumigatus: https://www.ncbi.nlm.nih.gov/sra/?term=SRX5000026
# * SRX6981195, SRX389628, grown in a murine chronic wound infection: https://www.ncbi.nlm.nih.gov/sra/?term=SRX69811951
#
# Looks like there are some samples that are isolated from wounds and others that are related to specific culturing or stress conditions like antibiotic treatment, temperature, co-culture with a fungus.
