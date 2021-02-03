
# coding: utf-8

# # Expression activity
# 
# This notebook examines the average expression activity of genes in order to:
# 1. Verify that accessory genes (i.e. PAO1 genes) mostly 0
# 2. Look at the data to determine if there are additional considerations to be made for our downstream analysis

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from textwrap import fill

from core_acc_modules import paths_corr, utils


# In[2]:


# Load gene expression data
expression_df = pd.read_csv(paths_corr.PAO1_GE, sep='\t', index_col=0, header=0)


# In[3]:


# Load metadata
pao1_metadata = pd.read_csv(paths_corr.PAO1_METADATA, index_col=0, header=0)
pa14_metadata = pd.read_csv(paths_corr.PA14_METADATA, index_col=0, header=0)


# ## Format gene ids

# In[4]:


# Replace sequencing gene ids with PA ids
pao1_fasta_file = paths_corr.PAO1_REF
seq_id_to_gene_id_pao1 = utils.dict_gene_num_to_ids(pao1_fasta_file)
expression_df.rename(mapper=seq_id_to_gene_id_pao1, axis="columns", inplace=True)

expression_df.head()


# ## Get core and accessory genes

# In[5]:


# Get mapping between PAO1 and PA14 genes using PAO1 reference
gene_annot_file = paths_corr.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(gene_annot_file, 'pao1')
gene_mapping_pao1.head()


# In[6]:


# Get mapping between PAO1 and PA14 genes using PA14 reference
gene_annot_file = paths_corr.GENE_PA14_ANNOT
gene_mapping_pa14 = utils.get_pao1_pa14_gene_map(gene_annot_file, 'pa14')
gene_mapping_pa14.head()


# In[7]:


# Get core genes: genes that have a homolog between PAO1 and PA14
core_pao1_genes, core_pa14_genes = utils.get_core_genes(gene_mapping_pao1,
                                                        gene_mapping_pa14,
                                                        False)
print(f"Number of PAO1 core genes: {len(core_pao1_genes)}")


# In[8]:


# Select only core genes that are included in my dataset
pao1_ref_genes = expression_df.columns
my_core_pao1_genes = list(set(core_pao1_genes).intersection(pao1_ref_genes))

print(f"Number of PAO1 core genes in my dataset: {len(my_core_pao1_genes)}")


# In[9]:


# Get PAO1-specific genes
pao1_acc = list(set(pao1_ref_genes) - set(my_core_pao1_genes))
print(f"Number of PAO1-specific genes: {len(pao1_acc)}")


# In[10]:


assert len(expression_df.columns) == len(my_core_pao1_genes) + len(pao1_acc)


# In[11]:


# Check that `get_pao1_pa14_gene_map` function is working as expected
assert("PA0053" not in core_pao1_genes and "PA0053" in pao1_acc)


# ## Get sample groups

# In[12]:


# Get PAO1 and PA14 sample ids
pao1_ids = list(pao1_metadata.index)
pa14_ids = list(pa14_metadata.index)


# ## Plot expression activity

# In[13]:


# Look at the expression of core genes in PAO1 and PA14 samples
expression_core_pao1 = expression_df.loc[pao1_ids, my_core_pao1_genes]
expression_core_pa14 = expression_df.loc[pa14_ids, my_core_pao1_genes]

expression_core_pao1_mean = expression_core_pao1.mean()
expression_core_pa14_mean = expression_core_pa14.mean()


# In[14]:


g = sns.distplot(expression_core_pao1_mean.values, kde=False)
g = sns.distplot(expression_core_pa14_mean.values, kde=False)
g.set_title('Histogram of mean gene expression for core genes')
g.set_xlabel('average TPM')
g.set_ylabel('Count')


# In[15]:


# Look at the expression of core genes in PAO1 and PA14 samples
expression_acc_pao1 = expression_df.loc[pao1_ids, pao1_acc]
expression_acc_pa14 = expression_df.loc[pa14_ids, pao1_acc]

expression_acc_pao1_mean = expression_acc_pao1.mean()
expression_acc_pa14_mean = expression_acc_pa14.mean()


# In[16]:


print(sum(expression_core_pao1_mean == 0))
print(sum(expression_acc_pao1_mean == 0))
print(sum(expression_core_pa14_mean == 0))
print(sum(expression_acc_pa14_mean == 0))


# In[17]:


f = sns.distplot(expression_acc_pao1_mean.values, kde=False)
f = sns.distplot(expression_acc_pa14_mean.values, kde=False)
f.set_title('Histogram of mean gene expression for PAO1-only genes')
f.set_xlabel('average TPM')
f.set_ylabel('Count')


# * Overall most genes have a TPM < 5000 and there is a very long tail in the distribution which is expected
# * Looks like there is a slight shift in PA14 samples have lower expression of PAO1-only genes, as expected. There is still nonzero expression for many PAO1-only genes in PA14 samples. How should we deal with this?
#   * Should we use presence/absence information to exclude certain samples from our correlation on a per-gene basis? This will be difficult because we will need to get this information somehow
#   * Should we only consider samples within PAO1 for these PAO1-only genes?
# * Looks like there is higher expression in core genes for both PAO1 and PA14 samples, as expected
# * There is a difference in the overall expression level of core vs PAO1-only genes. Since we plan to use Spearman correlation which is based on the ordering of values, I don't think the raw activity should effect the correlation.
