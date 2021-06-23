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

# # Module composition
#
# Given our co-acting gene modules, we will now calculate the composition of those modules - are modules predominantly core, predominantly accessory or mixed.
#
# The strategy we will use will be the following:
#
# 1. Given _N_ accessory genes, sample _N_ core genes 100 times
# 2. For each module in the network, compare the number of accessory genes in the module to the distribution of core gene counts in that module from your 100 samples
#     * if the number of accessory genes < 10th quantile of core genes then the module is mostly core
#     * if the number of accessory genes > 90th quantile of core genes then the module is mostly accessory
#     * else the module is mixed
#
#
# Explain fisher's exact test***
#
# Here we are comparing the ratio of core/acc within a given module to the ratio of core/acc outside that module.

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
from core_acc_modules import utils, paths

random.seed(1)
# -

# User params
method = "affinity"
# num_samples = 100

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

# ### Fisher's exact test

# Get list of modules ids
pao1_module_ids = pao1_membership["module id"].unique()
pa14_module_ids = pa14_membership["module id"].unique()

pao1_gene_group_composition = pd.DataFrame(
    index=pao1_module_ids, columns=["odds ratio", "p-value", "module label"]
)
pa14_gene_group_composition = pd.DataFrame(
    index=pa14_module_ids, columns=["odds ratio", "p-value", "module label"]
)

pao1_gene_group_composition.head()


def label_modules(
    module_ids_list, membership_df, core_genes_list, acc_genes_list, out_df
):
    all_genes = list(membership_df.index)

    for module_id in module_ids_list:
        # For each module create a table
        # -----|inside module|outside module
        # core | # genes     | # genes
        # acc  | # genes     | # genes

        # Find genes in module and outside module
        genes_in_module = list(
            membership_df[membership_df["module id"] == module_id].index
        )
        genes_outside_module = list(set(all_genes).difference(genes_in_module))

        # Get core and accessory in module
        core_genes_in_module = [
            gene for gene in genes_in_module if gene in core_genes_list
        ]
        acc_genes_in_module = [
            gene for gene in genes_in_module if gene in acc_genes_list
        ]

        # Get core and accessory genes outside of module
        core_genes_outside_module = [
            gene for gene in genes_outside_module if gene in core_genes_list
        ]
        acc_genes_outside_module = [
            gene for gene in genes_outside_module if gene in acc_genes_list
        ]

        # Check
        assert len(all_genes) == len(core_genes_in_module) + len(
            acc_genes_in_module
        ) + len(core_genes_outside_module) + len(acc_genes_outside_module)
        assert len(genes_in_module) == len(core_genes_in_module) + len(
            acc_genes_in_module
        )

        # Make contingency table
        observed_contingency_table = np.array(
            [
                [len(core_genes_in_module), len(core_genes_outside_module)],
                [len(acc_genes_in_module), len(acc_genes_outside_module)],
            ]
        )
        odds_ratio, pval = scipy.stats.fisher_exact(
            observed_contingency_table, alternative="two-sided"
        )

        # Fill in df
        out_df.loc[module_id, "odds ratio"] = odds_ratio
        out_df.loc[module_id, "p-value"] = pval

        # Do we want to include odds ratio >1 and significant p-value???
        if odds_ratio > 1 and pval < 0.05:
            out_df.loc[module_id, "module label"] = "mostly core"
        elif odds_ratio < 1 and pval < 0.05:
            out_df.loc[module_id, "module label"] = "mostly accessory"
        else:
            out_df.loc[module_id, "module label"] = "mixed"

    return out_df


# %%time
# Get labels for PAO1 compendium
pao1_module_labels = label_modules(
    pao1_module_ids, pao1_membership, pao1_core, pao1_acc, pao1_gene_group_composition
)

# %%time
# Get labels of PA14 compendium
pa14_module_labels = label_modules(
    pa14_module_ids, pa14_membership, pa14_core, pa14_acc, pa14_gene_group_composition
)

pao1_module_labels.head()

# ### Get matched number of core genes

"""# Sample same number of accessory genes from core genes `num_samples` times
# Store sets of core genes in a df
ls_core_pao1_samples = []
ls_core_pa14_samples = []

num_pao1_acc = len(pao1_acc)
num_pa14_acc = len(pa14_acc)

for i in range(num_samples):
    pao1_sample = random.sample(pao1_core, num_pao1_acc)
    pa14_sample = random.sample(pa14_core, num_pa14_acc)

    ls_core_pao1_samples.append(pao1_sample)
    ls_core_pa14_samples.append(pa14_sample)"""

"""assert len(ls_core_pao1_samples) == num_samples
assert len(ls_core_pa14_samples) == num_samples"""

"""assert len(ls_core_pao1_samples[0]) == num_pao1_acc
assert len(ls_core_pa14_samples[0]) == num_pa14_acc"""

# ### Calculate composition of modules

"""# Get list of modules ids
pao1_module_ids = pao1_membership["module id"].unique()
pa14_module_ids = pa14_membership["module id"].unique()"""

"""print(len(pao1_module_ids))
print(len(pa14_module_ids))"""

"""# For each module get the distribution of number of core genes
pao1_gene_group_composition = pd.DataFrame(
    index=pao1_module_ids, columns=range(0, num_samples)
)
pa14_gene_group_composition = pd.DataFrame(
    index=pa14_module_ids, columns=range(0, num_samples)
)


def get_module_composition(membership_df, core_samples_list, acc_genes, composition_df):
    # Get number of core genes from each sampling per module
    for i in range(len(core_samples_list)):
        num_core_genes = membership_df.loc[core_samples_list[i]][
            "module id"
        ].value_counts()
        composition_df[i] = num_core_genes

    # Get the number of accessory genes per module
    num_acc_genes = membership_df.loc[acc_genes]["module id"].value_counts()
    composition_df["acc"] = num_acc_genes

    composition_df = composition_df.fillna(0)

    return composition_df"""

"""pao1_gene_group_composition_processed = get_module_composition(
    pao1_membership, ls_core_pao1_samples, pao1_acc, pao1_gene_group_composition
)
pa14_gene_group_composition_processed = get_module_composition(
    pa14_membership, ls_core_pa14_samples, pa14_acc, pa14_gene_group_composition
)"""

# +
# pao1_gene_group_composition_processed.head()

# +
# pa14_gene_group_composition_processed.head()
# -

"""# For a given module,
# If the number of accessory genes < 10th quantile of core genes then the module is core
# If the number of accessory genes > 90th quantile of core genes then the module is accessory
# Else the module is mixed
def label_modules(module_composition_df):
    core_composition = module_composition_df.drop("acc", axis=1)

    module_composition_df["module label"] = "mixed"
    module_composition_df.loc[
        module_composition_df["acc"] < core_composition.quantile(0.1, axis=1),
        "module label",
    ] = "mostly core"
    module_composition_df.loc[
        module_composition_df["acc"] > core_composition.quantile(0.9, axis=1),
        "module label",
    ] = "mostly accessory"

    return module_composition_df"""

"""pao1_module_labels = label_modules(pao1_gene_group_composition_processed)
pa14_module_labels = label_modules(pa14_gene_group_composition_processed)"""

# +
# pao1_module_labels.head()

# +
# pa14_module_labels.head()
# -

# ### Examine module composition

pao1_module_labels["module label"].value_counts()

pa14_module_labels["module label"].value_counts()

# Add size of modules to df
pao1_module_labels["size"] = pao1_membership["module id"].value_counts()
pa14_module_labels["size"] = pa14_membership["module id"].value_counts()

pao1_module_labels.head()

pa14_module_labels.head()

# +
# Size distributions of PAO1 modules
f1 = sns.displot(
    pao1_module_labels.loc[pao1_module_labels["module label"] == "mixed", "size"]
)
plt.title("Size distribution of PAO1 mixed modules")

f2 = sns.displot(
    pao1_module_labels.loc[
        pao1_module_labels["module label"] == "mostly accessory", "size"
    ]
)
plt.title("Size distribution of PAO1 mostly accessory modules")

f3 = sns.displot(
    pao1_module_labels.loc[pao1_module_labels["module label"] == "mostly core", "size"]
)
plt.title("Size distribution of PAO1 mostly core modules")

# +
# Size distributions of PA14 moduels
g1 = sns.displot(
    pa14_module_labels.loc[pa14_module_labels["module label"] == "mixed", "size"]
)
plt.title("Size distribution of PA14 mixed modules")

g2 = sns.displot(
    pa14_module_labels.loc[
        pa14_module_labels["module label"] == "mostly accessory", "size"
    ]
)
plt.title("Size distribution of PA14 mostly accessory modules")

g3 = sns.displot(
    pa14_module_labels.loc[pa14_module_labels["module label"] == "mostly core", "size"]
)
plt.title("Size distribution of PA14 mostly core modules")
# -

# ### Lookup which genes are in which module
#
# TO DO
# * Move this into a different notebook?
# * Add annotations for KEGG, GO, gene name here

"""pao1_module_labels_truncated = pao1_module_labels[["module label", "size"]]
pa14_module_labels_truncated = pa14_module_labels[["module label", "size"]]"""

"""# Map genes to modules to labels
pao1_gene_module_labels = pao1_membership.merge(
    pao1_module_labels_truncated, left_on="module id", right_index=True
)
pa14_gene_module_labels = pa14_membership.merge(
    pa14_module_labels_truncated, left_on="module id", right_index=True
)"""

"""# Save
pao1_gene_module_labels.to_csv(
    os.path.join(
        paths.LOCAL_DATA_DIR, "pao1_gene_module_labels.tsv"
    ),
    sep="\t"
)
pa14_gene_module_labels.to_csv(
    os.path.join(
        paths.LOCAL_DATA_DIR, "pa14_gene_module_labels.tsv"
    ),
    sep="\t"
)"""

# +
# TO DO: Compare the composition across partitions
# Save the matrix that maps module id to module label per partition
# In another notebook look through module labels per partition
# Calculate consistency of labeling using ARI for each pair of partitions
# -

# **Takeaway:**
# * Most modules are mixed, some are mostly accessory. Only PA14 compendium have some mostly core modules
# * PAO1 mixed and mostly accessory modules have similar sizes (~10 genes)
# * PA14 mixed modules are ~10 genes, mostly accessory modules are smaller ~5 genes, mostly core modules are larger ~ 20 genes
