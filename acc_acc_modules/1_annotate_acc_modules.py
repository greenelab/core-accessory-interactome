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

# # Add annotations
#
# This notebook takes the dataframe with information about module composition and their labels and adds additional annotations including:
#
# 1. Which gene is contained within the modules (both gene id and gene name)
# 2. Baseline expression and expression in some context of interest
# 3. How clustered the module is on the genome
# 4. KEGG pathways that genes are found in
# 5. GO pathways genes are found in
# 6. Regulon/operon genes are found in
#
# All this information will help _P. aeruginosa_ experiments filter and determine which module might be interesting to explore.

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import scipy
import pandas as pd
import numpy as np
from itertools import product
from core_acc_modules import paths

random.seed(1)
# -

# User param
method = "affinity"

# +
# Import gene memberships
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_modules_{method}_acc.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_modules_{method}_acc.tsv"
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", index_col=0, header=0)
pa14_membership = pd.read_csv(pa14_membership_filename, sep="\t", index_col=0, header=0)

# +
# Import gene metadata
pao1_gene_annot_filename = paths.GENE_PAO1_ANNOT
pa14_gene_annot_filename = paths.GENE_PA14_ANNOT

pao1_gene_annot = pd.read_csv(pao1_gene_annot_filename, index_col=0, header=0)
pa14_gene_annot = pd.read_csv(pa14_gene_annot_filename, index_col=0, header=0)
# -

# Import metadata of samples
metadata_filename = paths.SAMPLE_METADATA

pao1_gene_annot = pao1_gene_annot["Name"].to_frame("gene name")
pa14_gene_annot = pa14_gene_annot["Name"].to_frame("gene name")

# ## Add gene names

# Add gene names
pao1_gene_module_labels = pao1_membership.merge(
    pao1_gene_annot, left_index=True, right_index=True
)
pa14_gene_module_labels = pa14_membership.merge(
    pa14_gene_annot, left_index=True, right_index=True
)

print(pao1_gene_module_labels.shape)
pao1_gene_module_labels.head()

print(pa14_gene_module_labels.shape)
pa14_gene_module_labels.head()

# ## Add expression information
#
# 1. What is the baseline level of expression for each gene in the module?
# 2. What is the expression level of genes in a clinical context (i.e. clinical samples)?

# Read in expression data
pao1_compendium = pd.read_csv(paths.PAO1_COMPENDIUM, sep="\t", index_col=0)
pa14_compendium = pd.read_csv(paths.PA14_COMPENDIUM, sep="\t", index_col=0)

pao1_compendium.head()

# Calculate median expression across all samples
pao1_median_all = pao1_compendium.median().to_frame("median expression")
pa14_median_all = pa14_compendium.median().to_frame("median expression")

pao1_median_all.head()


# Select subset of samples and calculate the median expression across that subset of samples
# TO DO: Move this into utils
def get_sample_ids(
    metadata_filename, experiment_colname, sample_colname, experiment_id
):
    """
    Returns sample ids (found in gene expression df) associated with
    a given list of experiment ids (found in the metadata)

    Arguments
    ----------
    metadata_filename: str
        File containing metadata
    experiment_colname: str
        Column header that contains experiment id that maps expression data
        and metadata
    sample_colname: str
        Column header that contains sample id that maps expression data
        and metadata
    experiment_id: str
        Selected experiment id to grab samples from

    """
    # Read in metadata
    metadata = pd.read_csv(metadata_filename, header=0)
    metadata.set_index(experiment_colname, inplace=True)

    selected_metadata = metadata.loc[experiment_id]
    sample_ids = list(selected_metadata[sample_colname])

    return sample_ids


# +
# TO DO: Have Deb or Georgia select a study
# Looks like we removed many of the clinical isolates from this compendium with our strain binning
# selected_sample_ids = get_sample_ids(
#   metadata_filename, experiment_colname="SRA_study", sample_colname="Experiment", experiment_id="SRP063289")

# +
# Subset compendium
# subset_pao1_compendium = pao1_compendium.loc[selected_sample_ids]
# subset_pa14_compendium = pa14_compendium.loc[selected_sample_ids]

# +
# print(subset_pao1_compendium.shape)
# print(subset_pa14_compendium.shape)

# +
# pao1_median_subset = subset_pao1_compendium.median().to_frame("median subset expression")
# pa14_median_subset = subset_pa14_compendium.median().to_frame("median subset expression")
# -

# Add median expression to gene ids
pao1_gene_annot = pao1_gene_module_labels.merge(
    pao1_median_all, left_index=True, right_index=True, how="left"
)
pa14_gene_annot = pa14_gene_module_labels.merge(
    pa14_median_all, left_index=True, right_index=True, how="left"
)

# Add median subset expression to gene ids
"""pao1_gene_annot = pao1_gene_annot.merge(
    pao1_median_subset, left_index=True, right_index=True, how="left"
)
pa14_gene_annot = pa14_gene_annot.merge(
    pa14_median_subset, left_index=True, right_index=True, how="left"
)"""

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

print(pa14_gene_annot.shape)
pa14_gene_annot.head()


# ## Genome location information
#
# How far are genes from other genes in the same module?

# +
# TO DO: Move into scripts
# For genes in the same module, calculate the pairwise distance from each other
# Calculate the median pairwise distance to represent how spread the module is
# across the genome
# Other metrics?


def get_intra_module_dist(annot_df, pa_prefix):
    rows = []
    for grp_name, grp_df in annot_df.groupby("module id"):
        # print("module", grp_name)

        # Trim off "PA" and convert number to integer
        ids = grp_df.index

        # Convert trailing id numbers to floats
        num_ids = [float(_id.split(pa_prefix)[1]) for _id in ids]

        abs_dist = []
        for gene1, gene2 in product(num_ids, num_ids):
            if gene1 != gene2:
                dist = abs(gene1 - gene2)
                # print(gene1, gene2, dist)
                abs_dist.append(dist)

        median_module_dist = np.median(abs_dist)
        min_dist = np.min(abs_dist)
        max_dist = np.max(abs_dist)

        for _id in ids:
            rows.append(
                {
                    "gene id": _id,
                    "median pairwise dist": median_module_dist,
                    "min pairwise dist": min_dist,
                    "max pairwise dist": max_dist,
                }
            )

    module_dist = pd.DataFrame(rows)
    module_dist = module_dist.set_index("gene id")

    return module_dist


# -

pao1_module_dist = get_intra_module_dist(pao1_gene_annot, pa_prefix="PA")
pa14_module_dist = get_intra_module_dist(pa14_gene_annot, pa_prefix="PA14_")

pao1_module_dist.head(10)

pa14_module_dist.head(10)

# Add module distance to gene names
pao1_gene_annot = pao1_gene_annot.merge(
    pao1_module_dist, left_index=True, right_index=True, how="left"
)
pa14_gene_annot = pa14_gene_annot.merge(
    pa14_module_dist, left_index=True, right_index=True, how="left"
)

# ## Add KEGG pathways
#
# For each pathway, what genes are contained in it. This information is only available for PAO1.

# +
pao1_pathway_filename = "https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt"

pao1_pathways = pd.read_csv(pao1_pathway_filename, sep="\t", index_col=0, header=None)
# -

pao1_pathways[2] = pao1_pathways[2].str.split(";").apply(set)
pao1_pathways.index = pao1_pathways.index.str.split(" - ").str[0]
pao1_pathways.head()

gene_to_pathways_df = pd.DataFrame(
    index=pao1_gene_module_labels.index, columns=list(pao1_pathways.index)
)

# %%time
for gene in gene_to_pathways_df.index:
    gene_to_pathways_df.loc[gene] = [
        gene in pao1_pathways.loc[pathway, 2] for pathway in pao1_pathways.index
    ]

gene_to_pathways_df.head()

# Add gene name to pathway information
pao1_gene_annot = pao1_gene_annot.merge(
    gene_to_pathways_df, left_index=True, right_index=True, how="left"
)

# ## Import and format operon

pao1_operon_filename = paths.PAO1_OPERON
pa14_operon_filename = paths.PA14_OPERON

pao1_operon = pd.read_csv(pao1_operon_filename, index_col=0, header=0)
pa14_operon = pd.read_csv(pa14_operon_filename, index_col=0, header=0)

pao1_operon = pao1_operon.set_index("locus_tag")
pa14_operon = pa14_operon.set_index("locus_tag")

print(pao1_operon.shape)
pao1_operon.head()

print(pa14_operon.shape)
pa14_operon.head()

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

# Only include columns for gene id and operon_name
pao1_operon = pao1_operon["operon_name"].to_frame()
pa14_operon = pa14_operon["operon_name"].to_frame()

pao1_operon.head()

# Add operons to pathway annotations for PAO1
pao1_gene_annot = pao1_gene_annot.merge(
    pao1_operon, left_index=True, right_index=True, how="left"
)

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

# For PA14 we only have operon annotations
pa14_gene_annot = pa14_gene_annot.merge(
    pa14_operon, left_index=True, right_index=True, how="left"
)

# ## Add regulon
#
# For each regulon, what genes are contained in it. This information is only available for PAO1

# +
pao1_regulon_filename = "https://raw.githubusercontent.com/greenelab/core-accessory-interactome/6635c0e357c0172c2cebd0368648030e0ee4beaf/data/metadata/regulons_format.csv"

pao1_regulons = pd.read_csv(pao1_regulon_filename, index_col=0, header=0)
# -

pao1_regulons["Genes"] = pao1_regulons["Genes"].str.split(";").apply(set)

gene_to_regulons_df = pd.DataFrame(
    index=pao1_gene_module_labels.index, columns=list(pao1_regulons.index)
)

# %%time
for gene in gene_to_regulons_df.index:
    gene_to_regulons_df.loc[gene] = [
        gene in pao1_regulons.loc[regulon, "Genes"] for regulon in pao1_regulons.index
    ]

# Add regulons to other annotations
pao1_gene_annot = pao1_gene_annot.merge(
    gene_to_regulons_df, left_index=True, right_index=True, how="left"
)

print(pao1_gene_annot.shape)
pao1_gene_annot.head()

print(pa14_gene_annot.shape)
pa14_gene_annot.head()

# Save
pao1_gene_annot.to_csv(f"pao1_acc_gene_module_annotated_{method}.tsv", sep="\t")
pa14_gene_annot.to_csv(f"pa14_acc_gene_module_annotated_{method}.tsv", sep="\t")

# These annotations will be used to help _P. aeruginosa_ experts, like our collaborators, to determine what accessory-accessory modules to focus on.
#
#
# Note: Since genes can be in multiple KEGG pathways and regulons, each pathway and regulon are separate columns. Whereas operons are a single column since genes can belong to only a single operon.
