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

# # Calculate conditional probability
#
# Now that we have labels for which module is mostly core, mostly accessory or mixed, we can ask our original question: How is the expression of different gene groups coordinated? Specifically we can ask: Are accessory genes more likely to be co-expressed with other accessory genes?
#
# To answer this we can calculate the following conditional probability:
# $$
# Pr(\text{gene x in an acc module | gene x is acc gene}) = \frac{Pr(\text{gene x in acc module}\cap \text{gene x is acc gene})}{Pr(\text{gene x is acc gene})}
# $$

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import numpy as np
from core_acc_modules import utils, paths

# User params
method = "affinity"

# +
# Import module memberships -- import annotated df
pao1_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pao1_gene_module_annotated_{method}.tsv"
)
pa14_membership_filename = os.path.join(
    paths.LOCAL_DATA_DIR, f"pa14_gene_module_annotated_{method}.tsv"
)

pao1_membership = pd.read_csv(pao1_membership_filename, sep="\t", index_col=0, header=0)
pa14_membership = pd.read_csv(pa14_membership_filename, sep="\t", index_col=0, header=0)
# -

pao1_membership.head()

pa14_membership.head()

# ## Get core/accessory annotations

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

# ## Calculate likelihood
# $Pr(\text{gene x in acc module}\cap \text{gene x is acc gene})$ is the number of accessory genes in accessory modules

# Get "mostly accessory" modules
pao1_subset = pao1_membership[["module id", "module label", "num acc in module"]]
pao1_subset = pao1_subset.set_index("module id")
pao1_subset = pao1_subset.drop_duplicates()

# Get "mostly accessory" modules
pa14_subset = pa14_membership[["module id", "module label", "num acc in module"]]
pa14_subset = pa14_subset.set_index("module id")
pa14_subset = pa14_subset.drop_duplicates()

num_acc_gene_in_acc_mod_pao1 = pao1_subset.loc[
    pao1_subset["module label"] == "mostly accessory", "num acc in module"
].sum()
num_acc_gene_in_acc_mod_pa14 = pa14_subset.loc[
    pa14_subset["module label"] == "mostly accessory", "num acc in module"
].sum()

lik_pao1_acc = num_acc_gene_in_acc_mod_pao1 / len(pao1_acc)
lik_pa14_acc = num_acc_gene_in_acc_mod_pa14 / len(pa14_acc)

num_acc_gene_in_acc_mod_pao1

num_acc_gene_in_acc_mod_pa14

print(lik_pao1_acc)
print(lik_pa14_acc)

# ### Caclulate prior distribution
# $Pr(\text{gene x is acc gene})$ is the number of accessory genes divided by the total number of genes

pr_pao1_acc = len(pao1_acc) / (len(pao1_core) + len(pao1_acc))
pr_pa14_acc = len(pa14_acc) / (len(pa14_core) + len(pa14_acc))

print(pr_pao1_acc)
print(pr_pa14_acc)

# ## Calculate conditional probability

# +
pr_acc2acc_pao1 = lik_pao1_acc / pr_pao1_acc
pr_acc2acc_pa14 = lik_pa14_acc / pr_pa14_acc

print(
    f"Probability of accessory gene being co-expressed with another accessory gene in PAO1 is {pr_acc2acc_pao1}"
)
print(
    f"Probability of accessory gene being co-expressed with another accessory gene in PA14 is {pr_acc2acc_pa14}"
)
