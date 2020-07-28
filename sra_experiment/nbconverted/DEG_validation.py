
# coding: utf-8

# # Differential expression validation
# This notebook performs a differential expression (DE) analysis comparing PAO1 samples vs PA14 samples. We can compare our results with those published in the literature as an additional step to validate that our RNA-seq processing are reasonable.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd
from core_acc_modules import utils, paths
from rpy2.robjects import pandas2ri
pandas2ri.activate()


# In[2]:


# Load gene expression using PAO1 reference
expression_data = pd.read_csv(paths.PAO1_GE, sep='\t', header=0, index_col=0)
print(expression_data.shape)
expression_data.head()


# ### Get core genes

# In[3]:


# Get mapping between PAO1 and PA14 genes using PAO1 reference
gene_annot_file = paths.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(gene_annot_file, 'pao1')
gene_annot_file = paths.GENE_PA14_ANNOT
gene_mapping_pa14 = utils.get_pao1_pa14_gene_map(gene_annot_file, 'pa14')

core_pao1_genes, core_pa14_genes = utils.get_core_genes(gene_mapping_pao1,
                                                        gene_mapping_pa14,
                                                        False)
print(f"Number of PAO1 core genes: {len(core_pao1_genes)}")
print(f"Number of PA14 core genes: {len(core_pa14_genes)}")

expression_data = expression_data.reindex(columns=core_pao1_genes)
print(expression_data.shape)
expression_data.head()


# In[4]:


gene_mapping_pa14.head()


# In[5]:


# Save file
expression_data.to_csv(paths.PAO1_GE_DE, sep='\t')


# ### Differential expression analysis

# In[ ]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Run one time\nif (!requireNamespace("BiocManager", quietly = TRUE))\n    install.packages("BiocManager")\nBiocManager::install("limma")')


# In[ ]:


get_ipython().run_cell_magic('R', '', "library('limma')")


# In[ ]:


# Load files
metadata_file = paths.SAMPLE_ANNOT
expression_data_file = paths.PAO1_GE_DE
out_file = paths.DE_STATS


# In[ ]:


get_ipython().run_cell_magic('R', '-i metadata_file -i expression_data_file -i out_file', "source('../core_acc_modules/DE_analysis.R')\n\nget_DE_stats(metadata_file,\n             expression_data_file,\n             out_file)")


# In[ ]:


# Read in DE stats file
DE_stats = pd.read_csv(paths.DE_STATS, sep='\t', header=0, index_col=0)
print(DE_stats.shape)
DE_stats.head()
# Get number of DEGs
# Get list of DEGs


# In[ ]:


# Compare out results with publication
#https://jb.asm.org/content/201/21/e00362-19 found ~ 2K DEGs between 2 strains where QS genes were DEGs
# What have publications found?

