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

# # Relationships using genome distance vs expression distance
#
# In our attempt to label modules as "mostly core", "mostly accessory" or "mixed". We found that most modules were "mixed" and some were "mostly accessory". We noticed that there were many modules that had only core genes, yet were not found to be signficanlty "mostly core" based on our Fisher's exact test due to the small size of the modules as well as the large imbalance in the number of core:accessory genes.
#
# These small modules, which are due to operons, is biologically sensible but hard for us to apply statistics. We want to try to tease apart the co-expression relationships that are due to locations (i.e. being in the same operon) versus other functional reasons.
#
# Our strategy is the following:
# * For each accessory gene, is the 1-NN/2-NN/3-NN core or accessory? Same for core genes
# * For each accessory gene, is the highest correlated/2nd-highest correlated/3rd highest correlated gene core or accessory? Same for core genes.
#
# Then we can compare the trends seen in both

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from core_acc_modules import utils, paths

random.seed(1)

# +
# User params
method = "affinity"
offset_to_bin = 10

use_operon = True
sum_increment_to_use = 2

# Output filename
pao1_figure_filename = (
    "PAO1_genome_expression_relationships_2window_operon_corrected.svg"
)
pa14_figure_filename = (
    "PA14_genome_expression_relationships_2window_operon_corrected.svg"
)
# -

# ### Import gene ids

# +
# Import correlation matrix to get gene ids
pao1_corr_filename = paths.PAO1_CORR_RAW
pa14_corr_filename = paths.PA14_CORR_RAW

pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pa14_corr = pd.read_csv(pa14_corr_filename, sep="\t", index_col=0, header=0)
# -

# Make a dataframe with gene ids
pao1_membership = pd.DataFrame(data=[], index=pao1_corr.index)
print(pao1_membership.shape)
pao1_membership.head()

pa14_membership = pd.DataFrame(data=[], index=pa14_corr.index)
print(pa14_membership.shape)
pa14_membership.head()

# ### Import and format operon data

pao1_operon_filename = paths.PAO1_OPERON
pa14_operon_filename = paths.PA14_OPERON

pao1_operon = pd.read_csv(pao1_operon_filename, index_col=0, header=0)
pa14_operon = pd.read_csv(pa14_operon_filename, index_col=0, header=0)

pao1_operon.head()

pao1_operon = pao1_operon.set_index("locus_tag")
pa14_operon = pa14_operon.set_index("locus_tag")

print(pao1_operon.shape)
pao1_operon.head()

# +
# There are 247 PAO1 genes with multiple annotations
# This operon df contains annotations from predicted operons based on DOOR database
# predictions which make up the majority of the operons) as well as some that
# are curated (i.e. PseudoCAP)
# There are some that have multiple PseudoCAP annotations too

# Here we will keep the last PseudoCAP annotations
# Note: Do we want to discard these annotations all together
# or will these need to be carefully curated to determine which to keep?
# We will use the curated annotation here
pao1_operon = pao1_operon[~pao1_operon.index.duplicated(keep="last")]
pa14_operon = pa14_operon[~pa14_operon.index.duplicated(keep="last")]
# -

pao1_operon.head()

# Only include columns for gene id and operon_name
pao1_operon = pao1_operon["operon_name"].to_frame()
pa14_operon = pa14_operon["operon_name"].to_frame()

print(pao1_operon.shape)
pao1_operon.head()

if use_operon:
    pao1_operon_expression_to_use = pao1_operon
    pa14_operon_expression_to_use = pa14_operon
else:
    pao1_operon_expression_to_use = None
    pa14_operon_expression_to_use = None

# ### Map core/accessory labels to genes

# +
# Read in expression data
pao1_expression_filename = paths.PAO1_COMPENDIUM
pa14_expression_filename = paths.PA14_COMPENDIUM

pao1_expression = pd.read_csv(pao1_expression_filename, sep="\t", index_col=0, header=0)
pa14_expression = pd.read_csv(pa14_expression_filename, sep="\t", index_col=0, header=0)

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

pao1_membership.loc[pao1_core, "core/acc"] = "core"
pao1_membership.loc[pao1_acc, "core/acc"] = "acc"

# pa14_acc_shared = set(pa14_acc).intersection(pa14_gene_module_labels.index)
pa14_membership.loc[pa14_core, "core/acc"] = "core"
pa14_membership.loc[pa14_acc, "core/acc"] = "acc"

# Drop "module id" column
pao1_arr = pao1_membership
pa14_arr = pa14_membership

# Make sure to sort by gene id
# NOTE PA14 gene ids don't increment by 1, but by 10 or 20 are we missing some genes?
pao1_arr = pao1_arr.sort_index()
pa14_arr = pa14_arr.sort_index()

print(pao1_arr.shape)
pao1_arr.head()

pao1_arr.tail()

print(pa14_arr.shape)
pa14_arr.head()

pa14_arr.tail()

# +
# Fill in index of operon_df to include all genes
all_pao1_gene_ids = pao1_arr.index
all_pa14_gene_ids = pa14_arr.index

# Get missing gene ids
missing_pao1_gene_ids = set(all_pao1_gene_ids).difference(pao1_operon.index)
missing_pa14_gene_ids = set(all_pa14_gene_ids).difference(pa14_operon.index)

# Make dataframe with missing gene ids with np.nan values for operon_name
missing_pao1_gene_df = pd.DataFrame(
    data=np.nan, index=list(missing_pao1_gene_ids), columns=["operon_name"]
)
missing_pa14_gene_df = pd.DataFrame(
    data=np.nan, index=list(missing_pa14_gene_ids), columns=["operon_name"]
)

pao1_operon_genome_dist = pao1_operon.append(missing_pao1_gene_df)
pa14_operon_genome_dist = pa14_operon.append(missing_pa14_gene_df)

pao1_operon_genome_dist = pao1_operon_genome_dist.loc[all_pao1_gene_ids]
pa14_operon_genome_dist = pa14_operon_genome_dist.loc[all_pa14_gene_ids]
# -

print(pao1_operon_genome_dist.shape)
pao1_operon_genome_dist.tail()

print(pa14_operon_genome_dist.shape)
pa14_operon_genome_dist.tail()

if use_operon:
    pao1_operon_genome_to_use = pao1_operon_genome_dist
    pa14_operon_genome_to_use = pa14_operon_genome_dist
else:
    pao1_operon_genome_to_use = None
    pa14_operon_genome_to_use = None


# ## Find relationships using genome distance

def get_relationship_in_genome_space(core_acc_df, offset_to_bin, operon_df=None):
    gene_type_start = ["acc", "core"]
    gene_type_compare = ["acc", "core"]

    core_acc_df_len = len(core_acc_df)
    offset_max = core_acc_df_len - 1

    core_acc_df_pad = np.pad(
        core_acc_df["core/acc"], offset_max, "constant", constant_values="NA"
    )

    if operon_df is not None:
        operon_df_pad = np.pad(
            operon_df["operon_name"], offset_max, "constant", constant_values="NA"
        )

        assert core_acc_df.shape == operon_df.shape

    rows = []

    for gene_start in gene_type_start:
        for gene_compare in gene_type_compare:
            for offset in range(1, offset_max + 1):
                # Print statements to understand what is happening
                # print("start left", core_acc_df_pad[offset_max : core_acc_df_len + offset_max])
                # print("compare left", core_acc_df_pad[
                #             offset_max - offset : core_acc_df_len + offset_max - offset
                #         ])
                #
                # print("start left", operon_df_pad[offset_max : core_acc_df_len + offset_max])
                # print("compare left", operon_df_pad[
                #             offset_max - offset : core_acc_df_len + offset_max - offset
                #         ])

                # Compare left nearest neighbors: Are they core or accessory?
                core_acc_left = (
                    core_acc_df_pad[offset_max : core_acc_df_len + offset_max]
                    == gene_start
                ) & (
                    core_acc_df_pad[
                        offset_max - offset : core_acc_df_len + offset_max - offset
                    ]
                    == gene_compare
                )

                # Compare right nearest neighbors: Are they core or accessory?
                core_acc_right = (
                    core_acc_df_pad[offset_max : core_acc_df_len + offset_max]
                    == gene_start
                ) & (
                    core_acc_df_pad[
                        offset_max + offset : core_acc_df_len + offset_max + offset
                    ]
                    == gene_compare
                )

                if operon_df is not None:
                    # Compare left operons: Are the genes in the same operon?
                    operon_left = (
                        operon_df_pad[offset_max : core_acc_df_len + offset_max]
                        == operon_df_pad[
                            offset_max - offset : core_acc_df_len + offset_max - offset
                        ]
                    )

                    # Compare right operons: Are the genes in the same operon?
                    operon_right = (
                        operon_df_pad[offset_max : core_acc_df_len + offset_max]
                        == operon_df_pad[
                            offset_max + offset : core_acc_df_len + offset_max + offset
                        ]
                    )

                    # Sum all comparisons
                    counts = (core_acc_left & ~operon_left).sum() + (
                        core_acc_right & ~operon_right
                    ).sum()
                else:
                    counts = core_acc_left.sum() + core_acc_right.sum()

                rows.append(
                    {
                        "gene start": gene_start,
                        "gene compare": gene_compare,
                        "offset": offset,
                        "total": counts,
                    }
                )

    genome_dist_counts = pd.DataFrame(rows)

    # Bin distances above offset_to_bin
    long_dist = (
        genome_dist_counts.query("offset>@offset_to_bin")
        .groupby(["gene compare", "gene start"])["total"]
        .mean()
        .reset_index()
    )
    long_dist["offset"] = f"{offset_to_bin}+"
    genome_dist_counts = genome_dist_counts.query("offset<=@offset_to_bin").append(
        long_dist, ignore_index=True
    )

    return genome_dist_counts


genome_dist_counts_pao1 = get_relationship_in_genome_space(
    pao1_arr, offset_to_bin, pao1_operon_genome_to_use
)
genome_dist_counts_pa14 = get_relationship_in_genome_space(
    pa14_arr, offset_to_bin, pa14_operon_genome_to_use
)

genome_dist_counts_pao1.head()

genome_dist_counts_pa14.head()

# ## Find relationships using expression distance

# Correlation matrix files
pao1_corr_filename = paths.PAO1_CORR_LOG_SPELL
pa14_corr_filename = paths.PA14_CORR_LOG_SPELL

# Load correlation data
pao1_corr = pd.read_csv(pao1_corr_filename, sep="\t", index_col=0, header=0)
pa14_corr = pd.read_csv(pa14_corr_filename, sep="\t", index_col=0, header=0)


def get_relationship_in_expression_space(
    corr_df,
    genes_to_consider,
    gene_mapping_df,
    offset_to_bin,
    operon_df=None,
    sum_increment=1,
):
    # Get subset of genes
    corr_subset = corr_df.loc[genes_to_consider]

    rows = []
    for gene in corr_subset.index:

        if operon_df is not None:
            # This subset needs to be reset each iteration
            # since we are dropping columns below
            corr_subset = corr_df.loc[genes_to_consider]

            # Check if gene is found in an operon
            if gene in operon_df.index:
                # Find operons containing 'gene'
                group_name = operon_df.loc[gene, "operon_name"]

                # Dictionary format: pao1_operon_dict[operon_name] = [list of genes]
                operon_dict = operon_df.groupby("operon_name").groups
                co_operonic_genes = list(operon_dict[group_name])

                # Remove columns corresponding to co-operonic genes from the corr_subset
                co_operonic_genes_to_remove = list(
                    set(corr_subset.columns).intersection(co_operonic_genes)
                )
                corr_subset = corr_subset.drop(columns=co_operonic_genes_to_remove)

        offset_max = corr_subset.shape[1]

        top_corr_genes = list(corr_subset.loc[gene].nlargest(offset_max).index[1:])
        top_gene_labels = list(gene_mapping_df.loc[top_corr_genes, "core/acc"].values)
        rows.append(top_gene_labels)

    expression_dist_counts = pd.DataFrame(rows)

    # Count types of relationships
    expression_dist_counts_acc = (expression_dist_counts == "acc").sum().to_frame("acc")

    expression_dist_counts = expression_dist_counts_acc.join(
        (expression_dist_counts == "core").sum().to_frame("core")
    )
    if sum_increment > 1:
        expression_dist_counts = (
            expression_dist_counts.rolling(sum_increment)
            .sum()
            .iloc[::-1]
            .shift(1)
            .sort_index()
        )

    # Format counts for plotting
    expression_dist_counts = expression_dist_counts.melt(
        var_name="gene type", value_name="total", ignore_index=False
    )
    expression_dist_counts = expression_dist_counts.rename_axis("offset").reset_index()
    expression_dist_counts["offset"] = expression_dist_counts["offset"] + 1

    # Average counts for weaker correlation relationships
    weak_corr = (
        expression_dist_counts.query("offset>@offset_to_bin")
        .groupby("gene type")["total"]
        .mean()
        .to_frame()
    )
    weak_corr = weak_corr.reset_index()
    weak_corr["offset"] = f"+{offset_to_bin}"

    expression_dist_counts = expression_dist_counts.query(
        "offset<=@offset_to_bin"
    ).append(weak_corr, ignore_index=True)

    # Add proportion - How should we calculate proportion
    # Of all 1-NN, %accessory, %core
    # Or, of all the accessory genes, %are 1-NN, of all core genes, %are 1-NN
    # total_counts = expression_dist_counts.groupby("offset")["total"].sum()[1]
    # expression_dist_counts["proportion"] = expression_dist_counts["total"]/total_counts

    return expression_dist_counts


# %%time
expression_dist_counts_pao1_acc = get_relationship_in_expression_space(
    pao1_corr,
    pao1_acc,
    pao1_arr,
    offset_to_bin,
    pao1_operon_expression_to_use,
    sum_increment_to_use,
)
expression_dist_counts_pao1_core = get_relationship_in_expression_space(
    pao1_corr,
    pao1_core,
    pao1_arr,
    offset_to_bin,
    pao1_operon_expression_to_use,
    sum_increment_to_use,
)

# %%time
expression_dist_counts_pa14_acc = get_relationship_in_expression_space(
    pa14_corr,
    pa14_acc,
    pa14_arr,
    offset_to_bin,
    pa14_operon_expression_to_use,
    sum_increment_to_use,
)
expression_dist_counts_pa14_core = get_relationship_in_expression_space(
    pa14_corr,
    pa14_core,
    pa14_arr,
    offset_to_bin,
    pa14_operon_expression_to_use,
    sum_increment_to_use,
)

expression_dist_counts_pao1_acc.head()

expression_dist_counts_pao1_core.head()

expression_dist_counts_pa14_acc.head()

expression_dist_counts_pa14_core.head()

# ### Plot

# +
# Plot PAO1 trends
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))

fig = sns.barplot(
    data=genome_dist_counts_pao1[genome_dist_counts_pao1["gene start"] == "acc"],
    x="offset",
    y="total",
    hue="gene compare",
    ax=axes[0][0],
    palette=sns.color_palette("Paired"),
)
fig.legend_.remove()
fig.set_title("Starting with accessory gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Offset in genome space")

fig = sns.barplot(
    data=genome_dist_counts_pao1[genome_dist_counts_pao1["gene start"] == "core"],
    x="offset",
    y="total",
    hue="gene compare",
    ax=axes[1][0],
    palette=sns.color_palette("Paired"),
)
fig.legend_.remove()
fig.set_title("Starting with core gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Offset in genome space")

fig = sns.barplot(
    data=expression_dist_counts_pao1_acc,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes[0][1],
    palette=sns.color_palette("Paired"),
)
fig.legend_.remove()
fig.set_title("Starting with accessory gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Rank correlation in expression space")

fig = sns.barplot(
    data=expression_dist_counts_pao1_core,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes[1][1],
    palette=sns.color_palette("Paired"),
)
fig.legend_.remove()
fig.set_title("Starting with core gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Rank correlation in expression space")


# Note: We are creating a single global legend that apply
# to all the facets of this figure. To do this using
# matplotlib, we need to be a little creative here
# and add the legend to a new location that is applied
# to the figure and then remove the legend from the facet.
plt.legend(bbox_to_anchor=(1.05, 1.15), loc=2, borderaxespad=0.0)

# +
# Plot PA14 trends
fig2, axes2 = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))

fig2 = sns.barplot(
    data=genome_dist_counts_pa14[genome_dist_counts_pa14["gene start"] == "acc"],
    x="offset",
    y="total",
    hue="gene compare",
    ax=axes2[0][0],
    palette=sns.color_palette("Paired"),
)
fig2.legend_.remove()
fig2.set_title("Starting with accessory gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Offset in genome space")

fig2 = sns.barplot(
    data=expression_dist_counts_pa14_acc,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes2[0][1],
    palette=sns.color_palette("Paired"),
)
fig2.legend_.remove()
fig2.set_title("Starting with accessory gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Rank correlation in expression space")

fig2 = sns.barplot(
    data=genome_dist_counts_pa14[genome_dist_counts_pa14["gene start"] == "core"],
    x="offset",
    y="total",
    hue="gene compare",
    ax=axes2[1][0],
    palette=sns.color_palette("Paired"),
)
fig2.legend_.remove()
fig2.set_title("Starting with core gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Offset in genome space")

fig2 = sns.barplot(
    data=expression_dist_counts_pa14_core,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes2[1][1],
    palette=sns.color_palette("Paired"),
)
fig2.legend_.remove()
fig2.set_title("Starting with core gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Rank correlation in expression space")

plt.legend(bbox_to_anchor=(1.05, 1.15), loc=2, borderaxespad=0.0)

# +
# Save figures using operons*
# Save figures not using operons*
# Save figure with rolling sum and operons
# Save figure with rolling sum not using operons
fig.figure.savefig(
    pao1_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

fig2.figure.savefig(
    pa14_figure_filename,
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)
# -

# **Takeaway:**
#
# In genome space:
# * The closest non co-operonic neighbor to an accessory gene is a core gene for PAO1, but is an accessory gene for PA14
# * Including co-operonic genes, accessory genes are clustered together on the genome (i.e. clustered with other accessory genes compared to core genes), which is known: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3008168/
# * Starting with a core gene, at any distance you find core genes including and not including co-operonic genes.
# * The x-axis indicates (left panel of plots) is the number of nearest neighbors (NN), determined based on the gene id that is sorted.
#
# In expression space:
# * Accessory genes are more likely to be highly co-expressed with other accessory genes, even accessory genes farther away (some coordination outside of location). This relationship is stronger in PA14 than PAO1 (i.e. accessory genes are more highly correlated with other accessory genes at farther distances in PA14). I wonder why this is.
# * Core genes are highly correlated with other core genes, again, this may be due to the fact that there are so many more core genes.
# * The x-axis(right panels of plots) is the number of correlated genes (i.e. 1=top most correlated gene, 2 =2nd most correlated gene)
#
# Note:
# * The distance calculation in genome space is bi-directional (i.e. it considers genes that are 1-NN to the left and right) so each starting gene has a count of 1-2. Whereas in expression space, we are only looking for the most correlated gene, so each starting gene has a count of 1. This explains why the range for the expression space plots are roughly half of the genome distance plots. To correct for this, there is an option to sum the counts across the 2 most highly correlated genes. The trends are consistent using an increment of 1 vs 2.
# * There is a drop off in the `10+` column for genome distance. This is because as we consider farther away genes, we lose the ability to consider in both directions and there are fewer genes at that far a distance. Should we adjust for this in our calculation?
#
# Some things to note about this analysis that may need to be updated:
# * There are some operons that have multiple annotations, which one should we choose? Should we drop these from the analysis? Should we curate these to determine which ones?
# * When we sorted the gene ids, we found that PAO1 incremented by 1 and PA14 incremented by 10 or 20, are we missing genes for PA14? How much will this change our genome dist analysis?