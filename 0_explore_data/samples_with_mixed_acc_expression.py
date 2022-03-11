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
import plotnine as pn
from scripts import paths, utils

# Threshold of expression value to determine if a sample has high accessory expression in both
threshold = 25

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

print(pao1_expression_prebin.shape)
print(pa14_expression_prebin.shape)

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
pao1_acc_compendium = pao1_expression_prebin[pao1_acc]
pao1_acc_compendium["median acc expression"] = pao1_acc_compendium.median(axis=1)
# -

print(pao1_acc_compendium.shape)
pao1_acc_compendium.head()

# +
# Create accessory df for PA14 compendium
# accessory gene ids | median accessory expression | strain label

# PA14-only genes in PA14 compendium
pa14_acc_compendium = pa14_expression_prebin[pa14_acc]
pa14_acc_compendium["median acc expression"] = pa14_acc_compendium.median(axis=1)
# -

print(pa14_acc_compendium.shape)
pa14_acc_compendium.head()

# +
# Merge acc compendia
pao1_pa14_acc_compendium = pao1_acc_compendium.merge(
    pa14_acc_compendium,
    left_index=True,
    right_index=True,
    suffixes=["_pao1", "_pa14"],
)

print(pao1_pa14_acc_compendium.shape)
pao1_pa14_acc_compendium.head()
# -

# Get samples with high PAO1 and PA14 expression
hybrid_samples_df = pao1_pa14_acc_compendium[
    (pao1_pa14_acc_compendium["median acc expression_pao1"] > threshold)
    & (pao1_pa14_acc_compendium["median acc expression_pa14"] > threshold)
]

# Add label for samples that are hybrid
pao1_pa14_acc_compendium["label"] = ""
pao1_pa14_acc_compendium.loc[hybrid_samples_df.index, "label"] = "hybrid"

# +
# Confirm with plot
# Plot
fig3 = pn.ggplot(
    pao1_pa14_acc_compendium,
    pn.aes(x="median acc expression_pao1", y="median acc expression_pa14"),
)
fig3 += pn.geom_point(
    pn.aes(color="label"),
    alpha=0.3,
)
# fig3 += pn.scale_color_manual(values=colors)
fig3 += pn.labs(
    x="median expression of PAO1-only genes",
    y="median expression of PA14-only genes",
    title="Accessory gene expression for all samples",
    width=10,
)
fig3 += pn.theme_bw()
fig3 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=14),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=16),
    axis_text=pn.element_text(family="sans-serif", size=12),
    axis_title=pn.element_text(family="sans-serif", size=14),
)
fig3 += pn.guides(
    colour=pn.guide_legend(
        title="SRA strain type", override_aes={"alpha": 1, "size": 3}
    )
)

print(fig3)

# +
# Plot log-scaled
fig4 = pn.ggplot(
    pao1_pa14_acc_compendium,
    pn.aes(x="median acc expression_pao1", y="median acc expression_pa14"),
)
fig4 += pn.scales.scale_x_log10()
fig4 += pn.scales.scale_y_log10()
fig4 += pn.geom_point(pn.aes(color="label"), alpha=0.4)
# fig4 += pn.scale_color_manual(values=colors)
fig4 += pn.labs(
    x="median expression of PAO1-only genes",
    y="median expression of PA14-only genes",
    title="log10 MR normalized estimated counts of accessory genes",
)
fig4 += pn.theme_bw()
fig4 += pn.theme(
    legend_title_align="center",
    plot_background=pn.element_rect(fill="white"),
    legend_key=pn.element_rect(fill="white", colour="white"),
    legend_title=pn.element_text(family="sans-serif", size=15),
    legend_text=pn.element_text(family="sans-serif", size=12),
    plot_title=pn.element_text(family="sans-serif", size=15),
    axis_text=pn.element_text(family="sans-serif", size=10),
    axis_title=pn.element_text(family="sans-serif", size=12),
)
fig4 += pn.guides(colour=pn.guide_legend(override_aes={"alpha": 1}))


print(fig4)
# -

hybrid_samples_df.head(10)

# Who are these samples?
# * [SRX6437734](https://www.ncbi.nlm.nih.gov/sra/?term=SRX6437734) and SRX6437785: From an expteriment comparing transcriptional patterns in biofilm vs planktonic growth. This sample is a biofilm sample.
# * [SRX6911057](https://www.ncbi.nlm.nih.gov/sra/?term=SRX6911057) a clinical isolate.
#
# There are not too many that had a high expression in both but 2 are from the same experiment studying biofilm growth and the 1 is a clinical isolate.

hybrid_samples_df.to_csv("Pa_hybrid_samples.tsv", sep="\t")
