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

# ## Validation of network modules
#
# This notebook performs a couple of analyses to validate the co-expression modules generated:
# 1. We examine the size of modules
# 2. We examine how co-operonic/co-regulonic genes are clustered into a few modules. A similar analysis can be found [here](spell_vs_counts_experiment/1a_compare_SPELL_vs_counts_correlation.ipynb) comparing within vs between edges for a given regulon/geneset.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import utils, paths

np.random.seed(1)

# +
# User params
# Params to examine module size
clustering_method_list = ["dbscan", "hierarchal", "affinity"]

# Params for regulon/operon coverage
# Clustering method to examine regulon/operon coverage
# This method needs to be one of the ones listed above in `clustering_method_list`
method_toexamine = "affinity"

# Remove modules of this size or greater for analysis looking at coverage of regulon/operons
module_size_threshold = 1000

# Seed to use to randomly sample a matched-sized set of genes
# to compare against regulon/operon composition
sample_seed = 1

# Gene subset
gene_subset = "acc"

# How was data processed
processed = "spell"
# -

# ## Examine size of modules
#
# This will serve as a quick check that we are using reasonable clustering params in [2_get_network_communities.ipynb](2_get_network_communities.ipynb)

for method_name in clustering_method_list:
    print(f"Modules using clustering method: {method_name}")
    pao1_membership_filename = os.path.join(
        paths.LOCAL_DATA_DIR,
        f"pao1_modules_{method_name}_{gene_subset}_{processed}.tsv",
    )
    pa14_membership_filename = os.path.join(
        paths.LOCAL_DATA_DIR,
        f"pa14_modules_{method_name}_{gene_subset}_{processed}.tsv",
    )

    pao1_membership = pd.read_csv(
        pao1_membership_filename, sep="\t", header=0, index_col=0
    )
    pa14_membership = pd.read_csv(
        pa14_membership_filename, sep="\t", header=0, index_col=0
    )
    # Note: Sort module ids by occurence for plotting later
    pao1_membership.sort_values(by="module id", ascending=False, inplace=True)
    pa14_membership.sort_values(by="module id", ascending=False, inplace=True)

    print(pao1_membership["module id"].value_counts())
    print(pa14_membership["module id"].value_counts())


# plotting function
def plot_dist_modules(clustering_method_list, gene_subset):

    # Set up the matplotlib figure
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))
    axes = axes.ravel()

    for i in range(len(clustering_method_list)):
        pao1_membership_filename = os.path.join(
            paths.LOCAL_DATA_DIR,
            f"pao1_modules_{clustering_method_list[i]}_{gene_subset}_{processed}.tsv",
        )
        pa14_membership_filename = os.path.join(
            paths.LOCAL_DATA_DIR,
            f"pa14_modules_{clustering_method_list[i]}_{gene_subset}_{processed}.tsv",
        )

        pao1_membership = pd.read_csv(
            pao1_membership_filename, sep="\t", header=0, index_col=0
        )
        pa14_membership = pd.read_csv(
            pa14_membership_filename, sep="\t", header=0, index_col=0
        )

        fig = (
            pao1_membership["module id"]
            .value_counts()
            .sort_values(ascending=False)
            .reset_index()["module id"]
            .plot(ax=axes[i])
        )
        fig = (
            pa14_membership["module id"]
            .value_counts()
            .sort_values(ascending=False)
            .reset_index()["module id"]
            .plot(ax=axes[i])
        )

        fig.set_title(
            f"Histogram of size of modules using {clustering_method_list[i]}",
            fontsize=12,
        )
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles, ["PAO1", "PA14"], loc="upper right")


# Plot distribution of modules per clustering method
plot_dist_modules(clustering_method_list, gene_subset)

# **Takeaway:**
# Our expectation on size of modules would be 2-50 genes. Most operons have fewer than 10 genes and most regulons have fewer than 100 genes. Some examples that demonstrate the size of co-expression networks can be found in papers using ADAGE signatures to define modules:
# * Figure 5 in [eADAGE paper](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-017-1905-4.pdf)
# * Figure 7 in [Harty et al. paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6531624/)
# * Figure 2 in [Doing et al. paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008783)
#
# What did we find? Which method follows our expectation?
# * Looks like there is one large modules using DBSCAN clustering
# * There are more even sized modules using hierarchal clustering and affinity propogation so we will probably use one of these 2 methods.

# ## Examine composition of modules
#
# This is a negative control. We expect that genes within the same operon or regulon will cluster together (i.e. be within the same module). To test this we will calculate the probability that a pair of genes will be from the same module, given that they are both from the same regulon or operon. We will calculate this probability for each (module, regulon/operon) combination.
#
# _Some definitions:_
#
# [Operons](https://en.wikipedia.org/wiki/Operon#:~:text=An%20operon%20is%20made%20up,transcription%20of%20the%20structural%20genes.) are a group of genes that share a promoter (DNA sequence that is recognized by RNA polymerase and enables transcription) and an operator (DNA sequence that repressor binds to and blocks RNA polymerase). Therefore these group of genes are transcribed or turned off together (so we would expect a very high correlation amongst these genes)
#
# [Regulons](https://en.wikipedia.org/wiki/Regulon) are a group of genes that are regulated by the same regulatory protein. A regulon can be composed of multiple operons.

# +
# Load PAO1 regulon file
pao1_regulon_filename = paths.PAO1_REGULON

# Load operon files
pa14_operon_filename = paths.PA14_OPERON
pao1_operon_filename = paths.PAO1_OPERON

# Load membership for specific clustering method
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR,
    f"pao1_modules_{method_toexamine}_{gene_subset}_{processed}.tsv",
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR,
    f"pa14_modules_{method_toexamine}_{gene_subset}_{processed}.tsv",
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", header=0, index_col=0)
pa14_membership = pd.read_csv(pa14_membership_filename, sep="\t", header=0, index_col=0)
# -

print(pao1_membership.shape)
pao1_membership.head()

print(pa14_membership.shape)
pa14_membership.head()

# According to Jake relationships tend to be more meaningful if the module is smaller (e.g. if an operon with 5 genes is contained in a module consisting of 10 total genes, this seems more biologically/functionally meaningful than an operon with 5 genes contained in a module consisting of 500 genes).
#
# To correct for the single or couple very large modules, we will remove them from the analysis

# +
# Get module ids that exceed size limit
module_todrop = (
    pao1_membership["module id"]
    .value_counts()[
        (pao1_membership["module id"].value_counts() > module_size_threshold)
    ]
    .index
)

print(module_todrop)

# +
# Get genes to drop
genes_todrop = pao1_membership[pao1_membership["module id"].isin(module_todrop)].index

# Drop genes
pao1_membership = pao1_membership.drop(genes_todrop)
# -

# ### Format operon/regulon files
#
# * Remove genes from operons/regulons that don't have membership information
# * Make random list of genes with matched size
# * There are many single gene operons, we will remove these for this analysis

# +
# Read file
pao1_operon = pd.read_csv(pao1_operon_filename, index_col=0, header=0)
pa14_operon = pd.read_csv(pa14_operon_filename, index_col=0, header=0)

pao1_regulon = pd.read_csv(pao1_regulon_filename, index_col=0, header=0)

print(pao1_operon.shape)
pao1_operon.head()
# -

print(pa14_operon.shape)
pa14_operon.head()

print(pao1_regulon.shape)
pao1_regulon.head()

# Convert "Genes" column from str to list for regulon dataset
pao1_regulon["Genes"] = pao1_regulon["Genes"].str.split(";")

# Check if genes within operon/regulon have membership information
# Only keep genes that are found in "pao1_membership"
pao1_regulon["Genes_processed"] = pao1_regulon["Genes"].apply(
    lambda list_genes: [
        gene_id for gene_id in list_genes if gene_id in pao1_membership.index
    ]
)

# Add size of the operons
pao1_operon["size"] = pao1_operon["operon_name"].map(
    pao1_operon.groupby("operon_name")["locus_tag"].count()
)
pa14_operon["size"] = pa14_operon["operon_name"].map(
    pa14_operon.groupby("operon_name")["locus_tag"].count()
)

pao1_operon.head()

pa14_operon.head()

pao1_operon_len = []
for grp_name, grp_df in pao1_operon.groupby("operon_name"):
    pao1_operon_len.append(grp_df.shape[0])

pa14_operon_len = []
for grp_name, grp_df in pa14_operon.groupby("operon_name"):
    pa14_operon_len.append(grp_df.shape[0])

# Update length based on filtered gene list ("Genes_processed" column)
pao1_regulon["size"] = pao1_regulon["Genes_processed"].str.len()

# If number genes in operon are 1 then remove
# Drop operons and regulons that have 0 genes due to no module filtering
pao1_operon = pao1_operon.drop(pao1_operon.query("size<=1").index)
pa14_operon = pa14_operon.drop(pa14_operon.query("size<=1").index)
pao1_regulon = pao1_regulon.drop(pao1_regulon.query("size<=1").index)

print(pao1_operon.shape)
pao1_operon.head()

print(pa14_operon.shape)
pa14_operon.head()

print(pao1_regulon.shape)
pao1_regulon.head()

# ### Get regulon information using PA14 ids
#
# Note that we can only do this mapping for core genes

pa14_regulon = pao1_regulon.copy()

# Get mapping between PAO1 and PA14 genes using PAO1 reference
gene_annot_file = paths.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(gene_annot_file, "pao1")
gene_mapping_pao1.head()

pa14_regulon["Genes_processed"] = pa14_regulon["Genes_processed"].apply(
    lambda pao1_gene_list: gene_mapping_pao1.loc[pao1_gene_list, "PA14_ID"].values
)

# Update length based on filtered gene list ("Genes_processed" column)
pa14_regulon["size"] = pa14_regulon["Genes_processed"].str.len()

# +
# If genes didn't map then drop operon/regulon
pa14_regulon = pa14_regulon.drop(
    pa14_regulon[
        pa14_regulon["Genes_processed"].apply(
            lambda gene_list: pd.isna(gene_list).any()
        )
    ].index
)

print(pa14_regulon.shape)
# -

# If number genes in operon are 1 then remove
# Drop operons and regulons that have 0 genes due to no module filtering
pa14_regulon = pa14_regulon.drop(pa14_regulon.query("size<=1").index)

print(pa14_regulon.shape)
pa14_regulon.head()

# +
# Quick look at distribution of size of regulons and operons
fig, axes = plt.subplots(ncols=4, nrows=1, figsize=(12, 5))

sns.distplot(
    pao1_operon_len,
    label="PAO1 operon size",
    color="red",
    kde=False,
    ax=axes[0],
)

sns.distplot(
    pa14_operon_len,
    label="PA14 operon size",
    color="red",
    kde=False,
    ax=axes[1],
)

sns.distplot(
    pao1_regulon["size"],
    label="PAO1 regulon size",
    color="blue",
    kde=False,
    ax=axes[2],
)

sns.distplot(
    pa14_regulon["size"],
    label="PA14 regulon size",
    color="blue",
    kde=False,
    ax=axes[3],
)

fig.suptitle(
    "Histogram of size of operons/regulons after filtering by membership",
    fontsize=12,
)
axes[0].set_title("PAO1 operon size")
axes[1].set_title("PA14 operon size")
axes[2].set_title("PAO1 regulon size")
axes[3].set_title("PA14 regulon size")

axes[2].set_xlabel("")
axes[3].set_xlabel("")

fig.text(0.5, 0.04, "Module size", ha="center")
axes[0].set_ylabel("count")


# -

# ### Calculate the probabilities
#
# What is the probability that gene x and y are in the same module given that they are both from the same regulon/operon?
#
# Given: regulon A and module B
# $$
# Pr(x,y \in B|x,y \in A) = \frac{Pr(x,y \in B \cap x,y \in A)}{Pr(x,y \in A)}
# $$

def coverage_of_genesets(module_df, genesets_df, geneset_type):

    total_genes = module_df.shape[0]
    rows = []

    for module_id, module_genes_df in module_df.groupby("module id"):

        # Pr(x,y in operon/regulon A)
        if geneset_type == "operon":
            # Dictionary of probabilities for a given module mapped to all operons
            operon_probs = {}
            for operon_id, operon_df in genesets_df.groupby("operon_name"):
                num_geneset = operon_df.shape[0]
                pr_denom = (num_geneset / total_genes) ** 2

                # Pr(x,y in module B | x,y in operon A)
                operon_df = operon_df.set_index("locus_tag")
                shared_genes = set(operon_df.index).intersection(module_genes_df.index)

                pr_joint = (len(shared_genes) / total_genes) ** 2
                pr_final = pr_joint / pr_denom

                operon_probs[operon_id] = pr_final

            # Save only the best matched operon-module based on the p-value
            (best_operon_id, best_prob) = max(operon_probs.items(), key=lambda k: k[1])

            rows.append(
                {
                    "module id": module_id,
                    "operon id": best_operon_id,
                    "pr(x,y in module|x,y in operon)": best_prob,
                }
            )

        else:
            # Dictionary of probabilities for a given module mapped to all operons
            regulon_probs = {}
            for regulon_id in genesets_df.index:
                num_geneset = genesets_df.loc[regulon_id, "size"]
                pr_denom = (num_geneset / total_genes) ** 2

                # Pr(x,y in module B | x,y in operon A)
                shared_genes = set(
                    genesets_df.loc[regulon_id, "Genes_processed"]
                ).intersection(module_genes_df.index)
                pr_joint = (len(shared_genes) / total_genes) ** 2
                pr_final = (pr_joint) / pr_denom

            regulon_probs[regulon_id] = pr_final

            # Save only the best matched operon-module based on the p-value
            (best_regulon_id, best_prob) = max(
                regulon_probs.items(), key=lambda k: k[1]
            )

            rows.append(
                {
                    "module id": module_id,
                    "regulon id": best_regulon_id,
                    "pr(x,y in module|x,y in regulon)": best_prob,
                }
            )
    out_df = pd.DataFrame(rows)
    if geneset_type == "operon":
        assert (out_df["pr(x,y in module|x,y in operon)"] > 1).sum() == 0
    else:
        assert (out_df["pr(x,y in module|x,y in regulon)"] > 1).sum() == 0
    return out_df


# %%time
pao1_operon_prob = coverage_of_genesets(pao1_membership, pao1_operon, "operon")
pao1_operon_prob.head()

# %%time
pa14_operon_prob = coverage_of_genesets(pa14_membership, pa14_operon, "operon")
pa14_operon_prob.head()

# %%time
pao1_regulon_prob = coverage_of_genesets(pao1_membership, pao1_regulon, "regulon")
pao1_regulon_prob.head()

# %%time
if gene_subset == "core":
    pa14_regulon_prob = coverage_of_genesets(pa14_membership, pa14_regulon, "regulon")
    pa14_regulon_prob.head()

# +
# As a baseline make a membership df mapping genes to a shuffled set of module ids
pao1_membership_shuffle = pao1_membership.copy()
pao1_membership_shuffle["module id"] = np.random.permutation(
    pao1_membership_shuffle["module id"].values
)

pa14_membership_shuffle = pa14_membership.copy()
pa14_membership_shuffle["module id"] = np.random.permutation(
    pa14_membership_shuffle["module id"].values
)
# -

# %%time
pao1_operon_shuffle_prob = coverage_of_genesets(
    pao1_membership_shuffle, pao1_operon, "operon"
)
pao1_operon_shuffle_prob.head()

# %%time
pa14_operon_shuffle_prob = coverage_of_genesets(
    pa14_membership_shuffle, pao1_operon, "operon"
)
pa14_operon_shuffle_prob.head()

# %%time
pao1_regulon_shuffle_prob = coverage_of_genesets(
    pao1_membership_shuffle, pao1_regulon, "regulon"
)
pao1_regulon_shuffle_prob.head()

# %%time
if gene_subset == "core":
    pa14_regulon_shuffle_prob = coverage_of_genesets(
        pa14_membership_shuffle, pa14_regulon, "regulon"
    )
    pa14_regulon_shuffle_prob.head()

# ## Plot distribution of probabilities
#
# Can we identify those operons, regulons that have high probability of being in the same module

# +
# Plot operon coverage
# Note: We are only plotting the probabilities greater than 0 since there were many operons that had
# a 0 probability, likely due to the small size of the operons
# All probabilities for PA14 shuffled data is 0 which is why the plot is blank
fig_operon, axes = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))
bins_shared = np.linspace(0, 1)

fig_operon = sns.histplot(
    pao1_operon_prob.loc[
        pao1_operon_prob["pr(x,y in module|x,y in operon)"] > 0,
        "pr(x,y in module|x,y in operon)",
    ],
    bins=bins_shared,
    ax=axes[0],
    label="true",
)
fig_operon = sns.histplot(
    pa14_operon_prob.loc[
        pa14_operon_prob["pr(x,y in module|x,y in operon)"] > 0,
        "pr(x,y in module|x,y in operon)",
    ],
    bins=bins_shared,
    ax=axes[1],
    label="true",
)
fig_operon = sns.histplot(
    pao1_operon_shuffle_prob.loc[
        pao1_operon_shuffle_prob["pr(x,y in module|x,y in operon)"] > 0,
        "pr(x,y in module|x,y in operon)",
    ],
    bins=bins_shared,
    color="grey",
    ax=axes[0],
    label="shuffle",
)
fig_operon = sns.histplot(
    pa14_operon_shuffle_prob.loc[
        pa14_operon_shuffle_prob["pr(x,y in module|x,y in operon)"] > 0,
        "pr(x,y in module|x,y in operon)",
    ],
    bins=bins_shared,
    color="grey",
    ax=axes[1],
    label="shuffle",
)

axes[0].set_title("PAO1 operon coverage")
axes[1].set_title("PA14 operon coverage")

legend = axes[0].legend()
legend = axes[1].legend()
# -

pao1_operon_prob[pao1_operon_prob["pr(x,y in module|x,y in operon)"] < 0.5]

# What are these operons that have a low prbabilities
# Why are some operons not likely to be within the same module?
# Based on the describe statistics, there doesn't appear to be a clear reasoning
# Overall, it is good that most operons have a high probability of being found in the same module
low_pao1_prob_operons = pao1_operon_prob[
    pao1_operon_prob["pr(x,y in module|x,y in operon)"] < 0.5
]["operon id"]
high_pao1_prob_operons = pao1_operon_prob[
    pao1_operon_prob["pr(x,y in module|x,y in operon)"] >= 0.5
]["operon id"]
pao1_operon_tmp = pao1_operon.set_index("operon_name")

pao1_operon_tmp.loc[low_pao1_prob_operons].describe()

pao1_operon_tmp.loc[high_pao1_prob_operons].describe()

low_pa14_prob_operons = pa14_operon_prob[
    pa14_operon_prob["pr(x,y in module|x,y in operon)"] < 0.5
]["operon id"]
high_pa14_prob_operons = pa14_operon_prob[
    pa14_operon_prob["pr(x,y in module|x,y in operon)"] >= 0.5
]["operon id"]
pa14_operon_tmp = pa14_operon.set_index("operon_name")

pa14_operon_tmp.loc[low_pa14_prob_operons].sort_values(by="size")

pa14_operon_tmp.loc[high_pa14_prob_operons].sort_values(by="size")

pao1_operon_prob.describe()

pao1_operon_shuffle_prob.describe()

pa14_operon_prob.describe()

pa14_operon_shuffle_prob.describe()

# +
# Plot regulon coverage
fig_regulon, axes = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))
bins_shared = np.linspace(0, 1)

fig_regulon = sns.histplot(
    pao1_regulon_prob["pr(x,y in module|x,y in regulon)"],
    bins=bins_shared,
    ax=axes[0],
    label="true",
)
if gene_subset == "core":
    fig_regulon = sns.histplot(
        pa14_regulon_prob["pr(x,y in module|x,y in regulon)"],
        bins=bins_shared,
        ax=axes[1],
    )
    fig_regulon = sns.histplot(
        pa14_regulon_shuffle_prob["pr(x,y in module|x,y in regulon)"],
        bins=bins_shared,
        color="grey",
        ax=axes[1],
    )
fig_regulon = sns.histplot(
    pao1_regulon_shuffle_prob["pr(x,y in module|x,y in regulon)"],
    bins=bins_shared,
    color="grey",
    ax=axes[0],
    label="shuffle",
)

axes[0].set_title("PAO1 regulon coverage")
if gene_subset == "core":
    axes[1].set_title("PA14 regulon coverage")

legend = axes[0].legend()
# -

pao1_regulon_prob.describe()

pao1_regulon_shuffle_prob.describe()

if gene_subset == "core":
    pa14_regulon_prob.describe()

if gene_subset == "core":
    pa14_regulon_shuffle_prob.describe()

# **Takeaway:**
# There is a higher probability that given pair of genes that are from the same operon, that they are also from the same module, compared to a randomly shuffled set of module assignments. Although there are some operons with low probabilties, overall genes in most operons have a high probability of being found in the same module.
#
# We don't see as drastic of a skewing for the regulons, though the mean using the true module labels is slightly higher compared to the shuffle module labels.
#
# Overall, this demonstrated that operons are well captured in our correlation matrix. However, a more effective way to assess this can be found [here](spell_vs_counts_experiment/1a_compare_SPELL_vs_counts_correlation.ipynb). Since the size of the regulons with respect to the non-regulon genes is so different, dividing by the total number of genes will drown out any signal we have.
