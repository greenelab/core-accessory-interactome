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

# +
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from core_acc_modules import utils, paths

random.seed(1)
# -

# User params
method = "affinity"
num_samples = 100

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

# ### Get matched number of core genes

# +
# Sample same number of accessory genes from core genes `num_samples` times
# Store sets of core genes in a df
ls_core_pao1_samples = []
ls_core_pa14_samples = []

num_pao1_acc = len(pao1_acc)
num_pa14_acc = len(pa14_acc)

for i in range(num_samples):
    pao1_sample = random.sample(pao1_core, num_pao1_acc)
    pa14_sample = random.sample(pa14_core, num_pa14_acc)

    ls_core_pao1_samples.append(pao1_sample)
    ls_core_pa14_samples.append(pa14_sample)
# -

assert len(ls_core_pao1_samples) == num_samples
assert len(ls_core_pa14_samples) == num_samples

assert len(ls_core_pao1_samples[0]) == num_pao1_acc
assert len(ls_core_pa14_samples[0]) == num_pa14_acc

# ### Calculate composition of modules

# Get list of modules ids
pao1_module_ids = pao1_membership["module id"].unique()
pa14_module_ids = pa14_membership["module id"].unique()

print(len(pao1_module_ids))
print(len(pa14_module_ids))

# +
# For each module get the distribution of number of core genes
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

    return composition_df


# -

pao1_gene_group_composition_processed = get_module_composition(
    pao1_membership, ls_core_pao1_samples, pao1_acc, pao1_gene_group_composition
)
pa14_gene_group_composition_processed = get_module_composition(
    pa14_membership, ls_core_pa14_samples, pa14_acc, pa14_gene_group_composition
)

pao1_gene_group_composition_processed.head()

pa14_gene_group_composition_processed.head()


# For a given module,
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

    return module_composition_df


pao1_module_labels = label_modules(pao1_gene_group_composition_processed)
pa14_module_labels = label_modules(pa14_gene_group_composition_processed)

pao1_module_labels.head()

pa14_module_labels.head()

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
