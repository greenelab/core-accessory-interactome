
# coding: utf-8

# # Explore highly co-expressed genes
# In the previous [notebook](2_explore_data.ipynb) we observed that using 39 samples with 201 PAO1-specific genes, that the correlation of accessory-accessory genes is higher compared to the correlation of core-core and core-accessory genes.
# 
# Based on this finding, we want to know: *What can explain this difference in correlation distribution?*
# 
# This notebook performs a follow-up analysis. In particular this notebook performs a deeper examination of the correlation structure per group (core-core, core-accessory, accessory-accessory) by looking at the trends of the nearest neighbors (i.e. highly correlated genes) of each gene.

# In[1]:


import pandas as pd
import os
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from functions import calculations

np.random.seed(123)


# In[2]:


# Input
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

base_intermediate_dir = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "tmp")

core_gene_ids_file = os.path.join(
    base_intermediate_dir,
    "core_gene_ids.pickle")

acc_gene_ids_file = os.path.join(
    base_intermediate_dir,
    "acc_gene_ids.pickle")

real_all_corr_file = os.path.join(
    base_intermediate_dir,
    "real_all_corr.pickle")

shuffled_all_corr_file = os.path.join(
    base_intermediate_dir,
    "shuffled_all_corr.pickle")

# Import Pseudomonas operon annotations from ADAGE repo
# Original source of data is from DOOR
# https://github.com/greenelab/adage/blob/master/Genome_organization/operon_3.txt
# Operons containing at least 3 genes
operon_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "annotations",
    "DOOR_operon_3.txt")


# # Examine highly co-expressed gene clusters
# For each core gene we will:
# 1. Extract the number of genes that are highly co-expressed with it
# 2. Determine the ratio of co-expressed genes that are core vs accessory
# 
# Repeat this for each accessory gene

# In[3]:


# Define threshold for highly co-expressed genes
coexpression_threshold = 0.75


# ## Co-expression patterns in real data

# In[4]:


# Get co-expression patterns using real expression data
real_core_df, real_acc_df = calculations.get_coexpression_stats(real_all_corr_file,
                                                                operon_file,
                                                                core_gene_ids_file,
                                                                acc_gene_ids_file,
                                                                coexpression_threshold)


# In[5]:


real_core_df.head()


# In[6]:


real_acc_df.head()


# In[7]:


# Print statistics about core genes
print('For a given CORE gene and using a threshold of {} to define co-expression: \n'.
     format(coexpression_threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median(real_core_df['num_coexpressed_genes'])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(real_core_df['percent_coexpressed_core'])*100,
             np.median(real_core_df['percent_coexpressed_acc'])*100))

# Calculate the percent of co-expressed genes that are non co-operonic
percent_non_cooperonic_coexpressed_genes = (real_core_df['num_non_cooperonic_coexpressed_genes']/real_core_df['num_coexpressed_genes']).fillna(0)

print('- There is a median of {}% co-expressed genes that are NOT in a shared operon'.
      format(np.median(percent_non_cooperonic_coexpressed_genes)*100))
print('- Of the non-coperonic co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(real_core_df['percent_non_cooperonic_coexpressed_core'])*100,
             np.median(real_core_df['percent_non_cooperonic_coexpressed_acc'])*100))


# In[8]:


# Print statistics about core genes
print('For a given ACCESSORY gene and using a threshold of {} to define co-expression: \n'.
     format(coexpression_threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median(real_acc_df['num_coexpressed_genes'])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(real_acc_df['percent_coexpressed_core'])*100,
             np.median(real_acc_df['percent_coexpressed_acc'])*100))

# Calculate the percent of co-expressed genes that are non co-operonic
percent_non_cooperonic_coexpressed_genes = (real_acc_df['num_non_cooperonic_coexpressed_genes']/real_acc_df['num_coexpressed_genes']).fillna(0)

print('- There is a median of {}% co-expressed genes that are NOT in a shared operon'.
      format(np.median(percent_non_cooperonic_coexpressed_genes)*100))
print('- Of the non-coperonic co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(real_acc_df['percent_non_cooperonic_coexpressed_core'])*100,
             np.median(real_acc_df['percent_non_cooperonic_coexpressed_acc'])*100))


# ## Co-expression patterns in shuffled data

# In[9]:


# Get co-expression patterns using shuffled expression data (control)
shuffled_core_df, shuffled_acc_df = calculations.get_coexpression_stats(shuffled_all_corr_file,
                                                                operon_file,
                                                                core_gene_ids_file,
                                                                acc_gene_ids_file,
                                                                coexpression_threshold)

shuffled_core_df.head()


# In[10]:


shuffled_acc_df.head()


# In[11]:


# Print statistics about core genes
print('For a given CORE gene and using a threshold of {} to define co-expression: \n'.
     format(coexpression_threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median(shuffled_core_df['num_coexpressed_genes'])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(shuffled_core_df['percent_coexpressed_core'])*100,
             np.median(shuffled_core_df['percent_coexpressed_acc'])*100))

# Calculate the percent of co-expressed genes that are non co-operonic
percent_non_cooperonic_coexpressed_genes = (shuffled_core_df['num_non_cooperonic_coexpressed_genes']/shuffled_core_df['num_coexpressed_genes']).fillna(0)

print('- There is a median of {}% co-expressed genes that are NOT in a shared operon'.
      format(np.median(percent_non_cooperonic_coexpressed_genes)*100))
print('- Of the non-coperonic co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(shuffled_core_df['percent_non_cooperonic_coexpressed_core'])*100,
             np.median(shuffled_core_df['percent_non_cooperonic_coexpressed_acc'])*100))


# In[12]:


# Print statistics about core genes
print('For a given ACCESSORY gene and using a threshold of {} to define co-expression: \n'.
     format(coexpression_threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median(shuffled_acc_df['num_coexpressed_genes'])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(shuffled_acc_df['percent_coexpressed_core'])*100,
             np.median(shuffled_acc_df['percent_coexpressed_acc'])*100))

# Calculate the percent of co-expressed genes that are non co-operonic
percent_non_cooperonic_coexpressed_genes = (shuffled_acc_df['num_non_cooperonic_coexpressed_genes']/shuffled_acc_df['num_coexpressed_genes']).fillna(0)

print('- There is a median of {}% co-expressed genes that are NOT in a shared operon'.
      format(np.median(percent_non_cooperonic_coexpressed_genes)*100))
print('- Of the non-coperonic co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median(shuffled_acc_df['percent_non_cooperonic_coexpressed_core'])*100,
             np.median(shuffled_acc_df['percent_non_cooperonic_coexpressed_acc'])*100))


# # Plot trends in co-expression data

# In[13]:


# Distribution of number of co-expressed genes
#sns.distplot(core_stats_df['num_coexpressed_genes'].tolist(), label='core', color='red')
#sns.distplot(acc_stats_df['num_coexpressed_genes'].tolist(), label='accessory', color='blue')

#plt.legend(prop={'size': 12})
#plt.title('Distribution of number of co-expressed genes')
#plt.xlabel('Number of co-expressed genes')
#plt.ylabel('Density')


# **Observation:**
# * Looks like both core and accessory genes do not tend to be co-expressed with many genes
# * Core genes are co-expressed with a median of 18 genes
# * Accessory genes are co-expressed with a median of 22 genes
# * Using 0.75 threshold to define co-expression
# * If using a more strict definition...

# In[14]:


# Distribution of number of non-operonic co-expressed genes
#sns.distplot(core_stats_df['num_coexpressed_genes'].tolist(), label='core', color='red')
#sns.distplot(acc_stats_df['num_coexpressed_genes'].tolist(), label='accessory', color='blue')

#plt.legend(prop={'size': 12})
#plt.title('Distribution of number of co-expressed genes')
#plt.xlabel('Number of co-expressed genes')
#plt.ylabel('Density')


# # Plot distribution of core, accessory co-expressed genes

# In[15]:


# Distribution plot for percent of core co-expressed genes
#sns.distplot(core_stats_df['percent_coexpressed_core'].tolist(), label='core', color='red')
#sns.distplot(acc_stats_df['percent_coexpressed_core'].tolist(), label='accessory', color='blue')

#plt.legend(prop={'size': 12})
#plt.title('Distribution of core co-expressed genes')
#plt.xlabel('Percent of core co-expressed genes')
#plt.ylabel('Density')


# In[16]:


# Distribution plot for percent of accessory co-expressed genes
#sns.distplot(core_stats_df['percent_coexpressed_acc'].tolist(), label='core', color='red')
#sns.distplot(acc_stats_df['percent_coexpressed_acc'].tolist(), label='accessory', color='blue')

#plt.legend(prop={'size': 12})
#plt.title('Distribution of accessory co-expressed genes')
#plt.xlabel('Percent of accessory co-expressed genes')
#plt.ylabel('Density')


# **Observation:**
# * Core genes tend to be co-expressed with only other core genes
# * Accessory genes tend to be co-expressed with some percent of core genes and accessory genes
