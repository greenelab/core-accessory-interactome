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

# # Common DEGs and core-accessory genes
#
# This notebook looks to see if the common DEGs, identified using [SOPHIE](https://github.com/greenelab/generic-expression-patterns/blob/master/pseudomonas_analysis/2_identify_generic_genes_pathways.ipynb) are mostly core or accessory genes.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import matplotlib
from matplotlib_venn import venn2
from core_acc_modules import utils, paths

# ### Get common DEGs

# +
# Read in SOPHIE identified common DEGs (PAO1 gene ids)
common_DEGs_filename = paths.COMMON_DEGS

common_DEGs = pd.read_csv(common_DEGs_filename, sep="\t", index_col=0, header=0)
# -

print(common_DEGs.shape)
common_DEGs.head()

common_DEGs = common_DEGs.set_index("gene id")

# ### Get core/accessory annotations

# +
# Read in expression data
pao1_expression_filename = paths.PAO1_COMPENDIUM
pa14_expression_filename = paths.PA14_COMPENDIUM

pao1_expression = pd.read_csv(pao1_expression_filename, sep="\t", index_col=0, header=0)
pa14_expression = pd.read_csv(pa14_expression_filename, sep="\t", index_col=0, header=0)
# -

# Note: Core and accessory annotations are from [BACTOME](https://academic.oup.com/nar/article/47/D1/D716/5112984). Not all core genes are measured by our expression dataset ("my dataset") we're using, so there may be a difference in "Number of PAO1 core genes" (core genes from BACTOME) and "Number of PAO1 core genes in my dataset" (core genes that are found in my expression dataset.

# +
pao1_annot_filename = paths.GENE_PAO1_ANNOT
pa14_annot_filename = paths.GENE_PA14_ANNOT

core_acc_dict = utils.get_my_core_acc_genes(
    pao1_annot_filename, pa14_annot_filename, pao1_expression, pa14_expression
)
# -

pao1_core = core_acc_dict["core_pao1"]
pa14_core = core_acc_dict["core_pa14"]
pao1_acc = core_acc_dict["acc_pao1"]
pa14_acc = core_acc_dict["acc_pa14"]

# ### Are common genes mostly core or accessory?

# Get shared gene ids
shared_acc = common_DEGs.index.intersection(pao1_acc)
shared_core = common_DEGs.index.intersection(pao1_core)

common_DEGs.loc[shared_acc, "gene group"] = "accessory"
common_DEGs.loc[shared_core, "gene group"] = "core"

common_DEGs["gene group"].value_counts()

# Add gene name
pao1_gene_annot = pd.read_csv(pao1_annot_filename, index_col=0, header=0)
pao1_gene_annot = pao1_gene_annot["Name"].to_frame("gene name")

common_DEGs_label = common_DEGs.merge(
    pao1_gene_annot, left_index=True, right_index=True
)

common_DEGs_label

# ### Venn diagram

common_DEGs_set = set(common_DEGs.index)
pao1_core_set = set(pao1_core)
pao1_acc_set = set(pao1_acc)

# +
core_common_venn = venn2(
    [common_DEGs_set, pao1_core_set], set_labels=("common DEGs", "Core genes")
)

core_common_venn.get_patch_by_id("11").set_color("purple")
core_common_venn.get_patch_by_id("11").set_edgecolor("none")
core_common_venn.get_patch_by_id("11").set_alpha(0.3)
core_common_venn.get_patch_by_id("01").set_color("blue")
core_common_venn.get_patch_by_id("01").set_edgecolor("none")
core_common_venn.get_patch_by_id("01").set_alpha(0.3)

matplotlib.pyplot.savefig(
    "common_core_venn.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# +
acc_common_venn = venn2(
    [common_DEGs_set, pao1_acc_set], set_labels=("common DEGs", "Accessory genes")
)

acc_common_venn.get_patch_by_id("11").set_color("purple")
acc_common_venn.get_patch_by_id("11").set_edgecolor("none")
acc_common_venn.get_patch_by_id("11").set_alpha(0.3)
acc_common_venn.get_patch_by_id("01").set_color("blue")
acc_common_venn.get_patch_by_id("01").set_edgecolor("none")
acc_common_venn.get_patch_by_id("01").set_alpha(0.3)

matplotlib.pyplot.savefig(
    "common_acc_venn.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# Save
common_DEGs_label.to_csv("common_DEGs_gene_group_labeled.tsv", sep="\t")

# **Takeaway:**
#
# Looks like most common DEGs are core, as expected. It is thought that these core genes encode essential functions shared by all strains and so it would make sense that these core genes are also those commonly DEGs.
