#!/usr/bin/env python
# coding: utf-8

# # Differential expression analysis
# 
# This notebook performs Differential Expression analysis using experiment, PRJNA283002, associated with [this publication](https://pubmed.ncbi.nlm.nih.gov/26078448/). Here they characterized the Anr regulon by comparing WT vs anr mutants.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
from rpy2.robjects import pandas2ri
pandas2ri.activate()

from core_acc_modules import paths_corr, utils, DE_helper


# In[2]:


# Load gene expression data
expression_df = pd.read_csv(paths_corr.PAO1_GE, sep='\t', index_col=0, header=0)


# ## Select expression data for experiment and replace gene ids

# In[3]:


# Select expression data associated with PRJNA283002 experiment
sample_metadata = pd.read_csv(paths_corr.DE_METADATA, sep='\t', index_col=0, header=0)
select_sample_ids = list(sample_metadata.index)

select_expression_df = expression_df.loc[select_sample_ids]


# In[4]:


# Replace gene sequencing ids with PAO1 ids to help us interpret our findings
pao1_fasta_file = paths_corr.PAO1_REF

seq_id_to_gene_id_pao1 = utils.dict_gene_num_to_ids(pao1_fasta_file)

select_expression_df.rename(mapper=seq_id_to_gene_id_pao1, axis="columns", inplace=True)

select_expression_df.head()


# In[5]:


# Save selected expression data
select_expression_df.to_csv(paths_corr.SELECT_GE, sep="\t")


# ## DE analysis

# In[6]:


# Process data for DESeq
DE_helper.process_samples_for_DESeq(
    paths_corr.SELECT_GE,
    paths_corr.DE_METADATA,
    paths_corr.SELECT_GE_PROCESSED,
)


# In[7]:


# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(paths_corr.DE_STATS_DIR, exist_ok=True)


# In[8]:


# Convert python path objects for use by R in the next cell
metadata_filename = str(paths_corr.DE_METADATA)
processed_expression_filename = str(paths_corr.SELECT_GE_PROCESSED)
repo_dir = str(paths_corr.PROJECT_DIR)
out_filename = str(paths_corr.DE_STATS_OUTPUT)


# In[9]:


get_ipython().run_cell_magic('R', '-i metadata_filename -i processed_expression_filename -i out_filename -i repo_dir', '\nsource(paste0(repo_dir, \'/core_acc_modules/DE_analysis.R\'))\n\n# File created: "<local_dir>/DE_stats/DE_stats_template_data_<project_id>_real.txt"\nget_DE_stats_DESeq(\n    metadata_filename,\n    processed_expression_filename,\n    out_filename\n)')


# ## Compare results with publication

# In[10]:


# Get top DEGs
# Compare against publication
DE_stats = pd.read_csv(
    paths_corr.DE_STATS_OUTPUT, 
    sep="\t", 
    header=0, 
    index_col=0
)

selected_DE_stats = DE_stats[(abs(DE_stats['log2FoldChange'])>1)]
print(selected_DE_stats.shape)
selected_DE_stats


# In[11]:


published_DEGs = [
    "PA1557",
    "PA3928",
    "PA2119",
    "PA3847",
    "PA0515",
    "PA0513",
    "PA0512",
    "PA0510",
    "PA0521",
    "PA0522",
    "PA0525",
    "PA0526",
    "PA2126",
    "PA2127",
    "PA2133",
]


# In[12]:


selected_DE_stats.loc[published_DEGs]


# In[13]:


input_DESeq_data = pd.read_csv(processed_expression_filename, sep="\t", index_col=0, header=0)

input_DESeq_data[published_DEGs]


# The differential expression results can be found in [Figure 1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4524035/) of the paper. Spot checking it looks like the genes have consistent direction of logFC.
# 
# Note:
# * NaN's occur if the samples are all 0. Need to check why PA3847 is NaN. Setting filtering to False doesn't get rid of NaNs
# 
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
