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
# -

# User params
method = "affinity"
offset_to_bin = 10

# ### Import module memberships

# +
# Import module memberships
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{method}.tsv"
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", index_col=0, header=0)
pa14_membership = pd.read_csv(pa14_membership_filename, sep="\t", index_col=0, header=0)
# -

pao1_membership.head()

pa14_membership.head()

# ### Import and format operon data

pao1_operon_filename = paths.PAO1_OPERON
pa14_operon_filename = paths.PA14_OPERON

pao1_operon = pd.read_csv(pao1_operon_filename, index_col=0, header=0)
pa14_operon = pd.read_csv(pa14_operon_filename, index_col=0, header=0)

pao1_operon = pao1_operon.set_index("locus_tag")
pa14_operon = pa14_operon.set_index("locus_tag")

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
pao1_arr = pao1_membership.drop("module id", axis=1)
pa14_arr = pa14_membership.drop("module id", axis=1)

pao1_arr.head()

pa14_arr.head()


# ## Find relationships using genome distance

def get_relationship_in_genome_space(core_acc_df, offset_to_bin):
    gene_type_start = ["acc", "core"]
    gene_type_compare = ["acc", "core"]

    core_acc_df_len = len(core_acc_df)
    offset_max = core_acc_df_len - 1

    core_acc_df_pad = np.pad(
        core_acc_df["core/acc"], offset_max, "constant", constant_values="NA"
    )

    rows = []

    for gene_start in gene_type_start:
        for gene_compare in gene_type_compare:
            for offset in range(1, offset_max + 1):
                # Print statements to understand what is happening
                # Here we are comparing our core/....
                # print(gene_start)
                # print(gene_compare)
                # print(pao1_arr_pad[offset_max:pao1_arr_len])
                # print(pao1_arr_pad[offset_max-offset:pao1_arr_len-offset])
                counts = (
                    (
                        core_acc_df_pad[offset_max : core_acc_df_len + offset_max]
                        == gene_start
                    )
                    & (
                        core_acc_df_pad[
                            offset_max - offset : core_acc_df_len + offset_max - offset
                        ]
                        == gene_compare
                    )
                ).sum() + (
                    (
                        core_acc_df_pad[offset_max : core_acc_df_len + offset_max]
                        == gene_start
                    )
                    & (
                        core_acc_df_pad[
                            offset_max + offset : core_acc_df_len + offset_max + offset
                        ]
                        == gene_compare
                    )
                ).sum()
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


genome_dist_counts_pao1 = get_relationship_in_genome_space(pao1_arr, offset_to_bin)
genome_dist_counts_pa14 = get_relationship_in_genome_space(pa14_arr, offset_to_bin)

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

                # Check if there are multiple annotations for this gene
                # This operon df contains annotations from predicted operons based on DOOR database
                # predictions which make up the majority of the operons) as well as some that
                # are curated (i.e. PseudoCAP)
                # We will use the curated annotation
                if type(group_name) != str:
                    group_name = (
                        operon_df.loc[gene]
                        .query("source_database=='PseudoCAP'")["operon_name"]
                        .values[0]
                    )

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
    pao1_corr, pao1_acc, pao1_arr, offset_to_bin, pao1_operon
)
expression_dist_counts_pao1_core = get_relationship_in_expression_space(
    pao1_corr, pao1_core, pao1_arr, offset_to_bin, pao1_operon
)

# %%time
expression_dist_counts_pa14_acc = get_relationship_in_expression_space(
    pa14_corr,
    pa14_acc,
    pa14_arr,
    offset_to_bin,
)
expression_dist_counts_pa14_core = get_relationship_in_expression_space(
    pa14_corr,
    pa14_core,
    pa14_arr,
    offset_to_bin,
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
)
fig.set_title("Starting with accessory gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Offset in genome space")

fig = sns.barplot(
    data=genome_dist_counts_pao1[genome_dist_counts_pao1["gene start"] == "core"],
    x="offset",
    y="total",
    hue="gene compare",
    ax=axes[1][0],
)
fig.set_title("Starting with core gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Offset in genome space")

fig = sns.barplot(
    data=expression_dist_counts_pao1_acc,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes[0][1],
)
fig.set_title("Starting with accessory gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Rank correlation in expression space")

fig = sns.barplot(
    data=expression_dist_counts_pao1_core,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes[1][1],
)
fig.set_title("Starting with core gene PAO1")
fig.set_ylabel("Number of genes")
fig.set_xlabel("Rank correlation in expression space")

# +
# Plot PA14 trends
fig2, axes2 = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))

fig2 = sns.barplot(
    data=genome_dist_counts_pa14[genome_dist_counts_pa14["gene start"] == "acc"],
    x="offset",
    y="total",
    hue="gene compare",
    ax=axes2[0][0],
)
fig2.set_title("Starting with accessory gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Offset in genome space")

fig2 = sns.barplot(
    data=expression_dist_counts_pa14_acc,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes2[0][1],
)
fig2.set_title("Starting with accessory gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Rank correlation in expression space")

fig2 = sns.barplot(
    data=genome_dist_counts_pa14[genome_dist_counts_pa14["gene start"] == "core"],
    x="offset",
    y="total",
    hue="gene compare",
    ax=axes2[1][0],
)
fig2.set_title("Starting with core gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Offset in genome space")

fig2 = sns.barplot(
    data=expression_dist_counts_pa14_core,
    x="offset",
    y="total",
    hue="gene type",
    ax=axes2[1][1],
)
fig2.set_title("Starting with core gene PA14")
fig2.set_ylabel("Number of genes")
fig2.set_xlabel("Rank correlation in expression space")
# -

# **Takeaway:**
#
# In genome space:
# * Accessory genes are clustered together on the genome (i.e. clustered with other accessory genes compared to core genes), which is known: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3008168/
# * Starting with a core gene, at any distance you find core genes because there are so many core genes
#
# In expression space:
# * Accessory genes appear to be more highly co-expressed with other accessory genes. But this might be due to the accessory genes clustering together in the genome (i.e. they are found in the same operon).
# * Core genes are highly correlated with other core genes, again, this may be due to the fact that there are so many more core genes.
