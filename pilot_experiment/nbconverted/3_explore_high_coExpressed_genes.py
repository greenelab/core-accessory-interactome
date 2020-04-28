
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


# In[3]:


# Read correlation data
core_gene_ids = pickle.load(open(core_gene_ids_file, "rb"))
acc_gene_ids = pickle.load(open(acc_gene_ids_file, "rb"))
real_all_corr = pickle.load(open(real_all_corr_file, "rb"))
shuffled_all_corr = pickle.load(open(shuffled_all_corr_file, "rb"))


# In[4]:


# Read operon data
# Manually had to set names to be the max size operon
operon_data = pd.read_csv(
    operon_file,
    header=None,
    sep='\t',
    names=range(15)
)

num_operons = operon_data.shape[0]
operon_data.head()


# In[5]:


# Get all gene ids
all_gene_ids = list(real_all_corr.index)


# # Examine highly co-expressed gene clusters
# For each core gene we will:
# 1. Extract the number of genes that are highly co-expressed with it
# 2. Determine the ratio of co-expressed genes that are core vs accessory
# 
# Repeat this for each accessory gene

# In[6]:


# Define threshold for highly co-expressed genes
threshold = 0.75


# In[7]:


# Apply threshold to identify which genes are co-expressed
real_all_coexpressed = real_all_corr>threshold
shuffled_all_coexpressed = shuffled_all_corr>threshold


# ## Get co-expressed genes using real data

# In[8]:


real_all_coexpressed.head(10)


# In[9]:


# Get upper triangle of correlation matrix
real_all_coexpressed_triu = pd.DataFrame(data=np.triu(real_all_coexpressed,1),
                                         index=real_all_coexpressed.index,
                                         columns=real_all_coexpressed.columns)

real_all_coexpressed_triu.head(10)


# ## Get number of co-expressed genes

# In[10]:


# Get total number of genes that are co-expressed per gene
num_coexpressed_genes = real_all_coexpressed_triu.sum(axis=1)

num_coexpressed_genes.head()


# ## Get list of co-expressed genes per gene id
# ## Get number of core and accessory co-expressed genes in list

# In[11]:


# Given the list of co-expressed genes
# we want to differentiate between those that are core and those that are accessory

name_coexpressed_genes = {}
num_coexpressed_core_genes = {}
num_coexpressed_acc_genes = {}

all_gene_ids = list(real_all_coexpressed_triu.index)

for gene_id in all_gene_ids:
    # Get row of correlation matrix
    # The values in the row corresponds to if there exists a gene is co-expressed with the gene_id
    coexpressed_gene_values = real_all_coexpressed_triu.loc[gene_id]
    
    # Check that our calculations are consistent
    assert(num_coexpressed_genes[gene_id] == sum(coexpressed_gene_values))
    
    if num_coexpressed_genes[gene_id] > 0:
        # Get list of co-expressed genes
        lst_coexpressed_genes = list(coexpressed_gene_values[coexpressed_gene_values].index)
        name_coexpressed_genes[gene_id] = lst_coexpressed_genes
        
        # Get the number of co-expressed genes that are core, accessory
        num_core_genes = len(set(lst_coexpressed_genes).intersection(core_gene_ids))
        num_acc_genes = len(set(lst_coexpressed_genes).intersection(acc_gene_ids))
        num_coexpressed_core_genes[gene_id] = num_core_genes
        num_coexpressed_acc_genes[gene_id] = num_acc_genes
        
    else:
        name_coexpressed_genes[gene_id] = []
        num_coexpressed_core_genes[gene_id] = 0
        num_coexpressed_acc_genes[gene_id] = 0


# In[12]:


# Calculate ratio of core:accessory genes in co-expressed gene sets
coexpressed_core_prop = {}
coexpressed_acc_prop = {}

for gene_id in all_gene_ids:
    num_core_genes = num_coexpressed_core_genes[gene_id]
    num_acc_genes = num_coexpressed_acc_genes[gene_id]
    if (num_core_genes == 0 & num_acc_genes == 0):
        coexpressed_core_prop[gene_id] = 0
        coexpressed_acc_prop[gene_id] = 0
    else:
        coexpressed_core_prop[gene_id] = num_core_genes/(num_core_genes + num_acc_genes)
        coexpressed_acc_prop[gene_id] = num_acc_genes/(num_core_genes + num_acc_genes)


# In[13]:


# Spot check that counts for what genes are core and accessroy are correct
#for i in ['PA0648', 'PA0980', 'PA1510', 'PA2037', 'PA2756', 'PA2794', 'PA3223', 'PA3793', 'PA4295', 'PA4451', 'PA4940', 'PA5086', 'PA5087']:
#    print(i in acc_gene_ids)


# ### Examine co-operonic co-expressed genes

# In[14]:


# Associate operons with reference gene
name_cooperonic_genes = {}

for gene_id in all_gene_ids:
    # Get operons containing reference gene_id   
    # Search for gene_id in each operon
    operon_search = operon_data.where(operon_data == gene_id).dropna(how='all').dropna(axis=1)
    
    if operon_search.empty:
        name_cooperonic_genes[gene_id] = []
    else:
        row_id = operon_search.index[0]
        name_cooperonic_genes[gene_id] = list(operon_data.loc[row_id].dropna())
        
name_cooperonic_genes    


# In[15]:


# Compare co-expressed gene set and co-operonic genes per reference gene id
num_non_cooperonic_coexpressed_genes = {}
num_non_cooperonic_coexpressed_core_genes = {}
num_non_cooperonic_coexpressed_acc_genes = {}

for gene_id in all_gene_ids:
    # Get co-operonic gene list
    cooperonic_genes = name_cooperonic_genes[gene_id]
    
    # Get co-expressed gene list
    coexpressed_genes = name_coexpressed_genes[gene_id]
    
    # Find non co-operonic genes
    # Find genes that DO NOT intersect between co-operonic genes and co-expressed genes
    cooperonic_coexpressed_genes = set(coexpressed_genes).intersection(cooperonic_genes)
    
    non_cooperonic_coexpressed_genes = set(coexpressed_genes) - cooperonic_coexpressed_genes
    
    # Get number of non-co-operonic genes
    num_non_cooperonic_coexpressed_genes[gene_id] = len(non_cooperonic_coexpressed_genes)
    
    if num_non_cooperonic_coexpressed_genes[gene_id] > 0:        
        # Get the number of non co-operonic co-expressed genes that are core, accessory
        num_core_genes = len(non_cooperonic_coexpressed_genes.intersection(core_gene_ids))
        num_acc_genes = len(non_cooperonic_coexpressed_genes.intersection(acc_gene_ids))
        num_non_cooperonic_coexpressed_core_genes[gene_id] = num_core_genes
        num_non_cooperonic_coexpressed_acc_genes[gene_id] = num_acc_genes
        
    else:
        num_non_cooperonic_coexpressed_core_genes[gene_id] = 0
        num_non_cooperonic_coexpressed_acc_genes[gene_id] = 0


# In[16]:


# Calculate ratio of core:accessory genes in co-expressed gene sets
non_cooperonic_coexpressed_core_prop = {}
non_cooperonic_coexpressed_acc_prop = {}

for gene_id in all_gene_ids:
    num_core_genes = num_non_cooperonic_coexpressed_core_genes[gene_id]
    num_acc_genes = num_non_cooperonic_coexpressed_acc_genes[gene_id]
    if (num_core_genes == 0 & num_acc_genes == 0):
        non_cooperonic_coexpressed_core_prop[gene_id] = 0
        non_cooperonic_coexpressed_acc_prop[gene_id] = 0
    else:
        non_cooperonic_coexpressed_core_prop[gene_id] = num_core_genes/(num_core_genes + num_acc_genes)
        non_cooperonic_coexpressed_acc_prop[gene_id] = num_acc_genes/(num_core_genes + num_acc_genes)


# # Summary statistics

# In[17]:


# Core gene stats
core_stats_df = pd.DataFrame(data={'ref_gene':core_gene_ids,
                                  'num_coexpressed_genes':num_coexpressed_genes[core_gene_ids],
                                   'num_coexpressed_core': [num_coexpressed_core_genes[k] for k in core_gene_ids],
                                   'num_coexpressed_acc': [num_coexpressed_acc_genes[k] for k in core_gene_ids],
                                   'percent_coexpressed_core': [coexpressed_core_prop[k] for k in core_gene_ids],
                                   'percent_coexpressed_acc': [coexpressed_acc_prop[k] for k in core_gene_ids],
                                   'num_non_cooperonic_coexpressed_genes':[num_non_cooperonic_coexpressed_genes[k] 
                                                                           for k in core_gene_ids],
                                   'num_non_cooperonic_coexpressed_core': [num_non_cooperonic_coexpressed_core_genes[k] 
                                                                           for k in core_gene_ids],
                                   'num_non_cooperonic_coexpressed_acc': [num_non_cooperonic_coexpressed_acc_genes[k] 
                                                                          for k in core_gene_ids],
                                   'percent_non_cooperonic_coexpressed_core': [non_cooperonic_coexpressed_core_prop[k] 
                                                                               for k in core_gene_ids],
                                   'percent_non_cooperonic_coexpressed_acc': [non_cooperonic_coexpressed_acc_prop[k] 
                                                                              for k in core_gene_ids]
                                  }
                            )
core_stats_df.head()


# In[18]:


# Accessory gene stats
acc_stats_df = pd.DataFrame(data={'ref_gene':acc_gene_ids,
                                  'num_coexpressed_genes':num_coexpressed_genes[acc_gene_ids],
                                  'num_coexpressed_core': [num_coexpressed_core_genes[a] for a in acc_gene_ids],
                                  'num_coexpressed_acc': [num_coexpressed_acc_genes[a] for a in acc_gene_ids],
                                  'percent_coexpressed_core': [coexpressed_core_prop[a] for a in acc_gene_ids],
                                  'percent_coexpressed_acc': [coexpressed_acc_prop[a] for a in acc_gene_ids],
                                  'num_non_cooperonic_coexpressed_genes':[num_non_cooperonic_coexpressed_genes[a] 
                                                                          for a in acc_gene_ids],
                                  'num_non_cooperonic_coexpressed_core': [num_non_cooperonic_coexpressed_core_genes[a] 
                                                                          for a in acc_gene_ids],
                                  'num_non_cooperonic_coexpressed_acc': [num_non_cooperonic_coexpressed_acc_genes[a] 
                                                                         for a in acc_gene_ids],
                                  'percent_non_cooperonic_coexpressed_core': [non_cooperonic_coexpressed_core_prop[a] 
                                                                               for a in acc_gene_ids],
                                  'percent_non_cooperonic_coexpressed_acc': [non_cooperonic_coexpressed_acc_prop[a] 
                                                                              for a in acc_gene_ids]
                                  }
                            )
acc_stats_df.head()


# ### Core gene statistics

# In[19]:


# Print statistics about core genes
print('For a given CORE gene and using a threshold of {} to define co-expression: \n'.
     format(threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median([num_coexpressed_genes[k] for k in core_gene_ids])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median([coexpressed_core_prop[k] for k in core_gene_ids])*100,
             np.median([coexpressed_acc_prop[k] for k in core_gene_ids])*100))


# ### Accessory gene statistics

# In[20]:


# Print statistics about core genes
print('For a given ACCESSORY gene and using a threshold of {} to define co-expression: \n'.
     format(threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median([num_coexpressed_genes[a] for a in acc_gene_ids])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median([coexpressed_core_prop[a] for a in acc_gene_ids])*100,
             np.median([coexpressed_acc_prop[a] for a in acc_gene_ids])*100))


# ## Using shuffled data

# In[21]:


# Apply threshold to identify which genes are co-expressed
shuffled_all_coexpressed = shuffled_all_corr>threshold


# In[22]:


# Get upper triangle of correlation matrix
shuffled_all_coexpressed_triu = pd.DataFrame(data=np.triu(shuffled_all_coexpressed,1),
                                         index=shuffled_all_coexpressed.index,
                                         columns=shuffled_all_coexpressed.columns)

shuffled_all_coexpressed_triu.head(10)


# In[23]:


# Get total number of genes that are co-expressed per gene
num_coexpressed_genes = shuffled_all_coexpressed_triu.sum(axis=1)


# In[24]:


# Given the list of co-expressed genes
# we want to differentiate between those that are core and those that are accessory

name_coexpressed_genes = {}
num_coexpressed_core_genes = {}
num_coexpressed_acc_genes = {}

all_gene_ids = list(real_all_coexpressed_triu.index)

for gene_id in all_gene_ids:
    # Get row of correlation matrix
    # The values in the row corresponds to if there exists a gene is co-expressed with the gene_id
    coexpressed_gene_values = shuffled_all_coexpressed_triu.loc[gene_id]
    
    # Check that our calculations are consistent
    assert(num_coexpressed_genes[gene_id] == sum(coexpressed_gene_values))
    
    if num_coexpressed_genes[gene_id] > 0:
        # Get list of co-expressed genes
        lst_coexpressed_genes = list(coexpressed_gene_values[coexpressed_gene_values].index)
        name_coexpressed_genes[gene_id] = lst_coexpressed_genes
        
        # Get the number of co-expressed genes that are core, accessory
        num_core_genes = len(set(lst_coexpressed_genes).intersection(core_gene_ids))
        num_acc_genes = len(set(lst_coexpressed_genes).intersection(acc_gene_ids))
        num_coexpressed_core_genes[gene_id] = num_core_genes
        num_coexpressed_acc_genes[gene_id] = num_acc_genes
        
    else:
        name_coexpressed_genes[gene_id] = []
        num_coexpressed_core_genes[gene_id] = 0
        num_coexpressed_acc_genes[gene_id] = 0


# In[25]:


# Calculate ratio of core:accessory genes in co-expressed gene sets
coexpressed_core_prop = {}
coexpressed_acc_prop = {}

for gene_id in all_gene_ids:
    num_core_genes = num_coexpressed_core_genes[gene_id]
    num_acc_genes = num_coexpressed_acc_genes[gene_id]
    if (num_core_genes == 0 & num_acc_genes == 0):
        coexpressed_core_prop[gene_id] = 0
        coexpressed_acc_prop[gene_id] = 0
    else:
        coexpressed_core_prop[gene_id] = num_core_genes/(num_core_genes + num_acc_genes)
        coexpressed_acc_prop[gene_id] = num_acc_genes/(num_core_genes + num_acc_genes)


# In[27]:


# Compare co-expressed gene set and co-operonic genes per reference gene id
num_non_cooperonic_coexpressed_genes = {}
num_non_cooperonic_coexpressed_core_genes = {}
num_non_cooperonic_coexpressed_acc_genes = {}

for gene_id in all_gene_ids:
    # Get co-operonic gene list
    cooperonic_genes = name_cooperonic_genes[gene_id]
    
    # Get co-expressed gene list
    coexpressed_genes = name_coexpressed_genes[gene_id]
    
    # Find non co-operonic genes
    # Find genes that DO NOT intersect between co-operonic genes and co-expressed genes
    cooperonic_coexpressed_genes = set(coexpressed_genes).intersection(cooperonic_genes)
    
    non_cooperonic_coexpressed_genes = set(coexpressed_genes) - cooperonic_coexpressed_genes
    
    # Get number of non-co-operonic genes
    num_non_cooperonic_coexpressed_genes[gene_id] = len(non_cooperonic_coexpressed_genes)
    
    if num_non_cooperonic_coexpressed_genes[gene_id] > 0:        
        # Get the number of non co-operonic co-expressed genes that are core, accessory
        num_core_genes = len(non_cooperonic_coexpressed_genes.intersection(core_gene_ids))
        num_acc_genes = len(non_cooperonic_coexpressed_genes.intersection(acc_gene_ids))
        num_non_cooperonic_coexpressed_core_genes[gene_id] = num_core_genes
        num_non_cooperonic_coexpressed_acc_genes[gene_id] = num_acc_genes
        
    else:
        num_non_cooperonic_coexpressed_core_genes[gene_id] = 0
        num_non_cooperonic_coexpressed_acc_genes[gene_id] = 0


# In[28]:


# Calculate ratio of core:accessory genes in co-expressed gene sets
non_cooperonic_coexpressed_core_prop = {}
non_cooperonic_coexpressed_acc_prop = {}

for gene_id in all_gene_ids:
    num_core_genes = num_non_cooperonic_coexpressed_core_genes[gene_id]
    num_acc_genes = num_non_cooperonic_coexpressed_acc_genes[gene_id]
    if (num_core_genes == 0 & num_acc_genes == 0):
        non_cooperonic_coexpressed_core_prop[gene_id] = 0
        non_cooperonic_coexpressed_acc_prop[gene_id] = 0
    else:
        non_cooperonic_coexpressed_core_prop[gene_id] = num_core_genes/(num_core_genes + num_acc_genes)
        non_cooperonic_coexpressed_acc_prop[gene_id] = num_acc_genes/(num_core_genes + num_acc_genes)


# In[29]:


# Core gene stats
core_stats_df = pd.DataFrame(data={'ref_gene':core_gene_ids,
                                  'num_coexpressed_genes':num_coexpressed_genes[core_gene_ids],
                                   'num_coexpressed_core': [num_coexpressed_core_genes[k] for k in core_gene_ids],
                                   'num_coexpressed_acc': [num_coexpressed_acc_genes[k] for k in core_gene_ids],
                                   'percent_coexpressed_core': [coexpressed_core_prop[k] for k in core_gene_ids],
                                   'percent_coexpressed_acc': [coexpressed_acc_prop[k] for k in core_gene_ids],
                                   'num_non_cooperonic_coexpressed_genes':[num_non_cooperonic_coexpressed_genes[k] 
                                                                           for k in core_gene_ids],
                                   'num_non_cooperonic_coexpressed_core': [num_non_cooperonic_coexpressed_core_genes[k] 
                                                                           for k in core_gene_ids],
                                   'num_non_cooperonic_coexpressed_acc': [num_non_cooperonic_coexpressed_acc_genes[k] 
                                                                          for k in core_gene_ids],
                                   'percent_non_cooperonic_coexpressed_core': [non_cooperonic_coexpressed_core_prop[k] 
                                                                               for k in core_gene_ids],
                                   'percent_non_cooperonic_coexpressed_acc': [non_cooperonic_coexpressed_acc_prop[k] 
                                                                              for k in core_gene_ids]
                                  }
                            )
core_stats_df.head()


# In[30]:


# Accessory gene stats
acc_stats_df = pd.DataFrame(data={'ref_gene':acc_gene_ids,
                                  'num_coexpressed_genes':num_coexpressed_genes[acc_gene_ids],
                                  'num_coexpressed_core': [num_coexpressed_core_genes[a] for a in acc_gene_ids],
                                  'num_coexpressed_acc': [num_coexpressed_acc_genes[a] for a in acc_gene_ids],
                                  'percent_coexpressed_core': [coexpressed_core_prop[a] for a in acc_gene_ids],
                                  'percent_coexpressed_acc': [coexpressed_acc_prop[a] for a in acc_gene_ids],
                                  'num_non_cooperonic_coexpressed_genes':[num_non_cooperonic_coexpressed_genes[a] 
                                                                          for a in acc_gene_ids],
                                  'num_non_cooperonic_coexpressed_core': [num_non_cooperonic_coexpressed_core_genes[a] 
                                                                          for a in acc_gene_ids],
                                  'num_non_cooperonic_coexpressed_acc': [num_non_cooperonic_coexpressed_acc_genes[a] 
                                                                         for a in acc_gene_ids],
                                  'percent_non_cooperonic_coexpressed_core': [non_cooperonic_coexpressed_core_prop[a] 
                                                                               for a in acc_gene_ids],
                                  'percent_non_cooperonic_coexpressed_acc': [non_cooperonic_coexpressed_acc_prop[a] 
                                                                              for a in acc_gene_ids]
                                  }
                            )
acc_stats_df.head()


# In[31]:


# Print statistics about core genes
print('For a given CORE gene and using a threshold of {} to define co-expression: \n'.
     format(threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median([num_coexpressed_genes[k] for k in core_gene_ids])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median([coexpressed_core_prop[k] for k in core_gene_ids])*100,
             np.median([coexpressed_acc_prop[k] for k in core_gene_ids])*100))


# In[33]:


# Print statistics about core genes
print('For a given ACCESSORY gene and using a threshold of {} to define co-expression: \n'.
     format(threshold))
print('- There is a median of {} co-expressed  genes'.
      format(np.median([num_coexpressed_genes[a] for a in acc_gene_ids])))
print('- Of the co-expressed genes, the median percent of core genes is {}% and accessory genes is {}%'.
      format(np.median([coexpressed_core_prop[a] for a in acc_gene_ids])*100,
             np.median([coexpressed_acc_prop[a] for a in acc_gene_ids])*100))

print('- There is a median of {} co-expressed genes that are NOT in a shared operon')
# Print same for core, acc


# # Plot number of co-expressed genes

# In[ ]:


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

# In[ ]:


# Distribution of number of non-operonic co-expressed genes
#sns.distplot(core_stats_df['num_coexpressed_genes'].tolist(), label='core', color='red')
#sns.distplot(acc_stats_df['num_coexpressed_genes'].tolist(), label='accessory', color='blue')

#plt.legend(prop={'size': 12})
#plt.title('Distribution of number of co-expressed genes')
#plt.xlabel('Number of co-expressed genes')
#plt.ylabel('Density')


# # Plot distribution of core, accessory co-expressed genes

# In[ ]:


# Distribution plot for percent of core co-expressed genes
#sns.distplot(core_stats_df['percent_coexpressed_core'].tolist(), label='core', color='red')
#sns.distplot(acc_stats_df['percent_coexpressed_core'].tolist(), label='accessory', color='blue')

#plt.legend(prop={'size': 12})
#plt.title('Distribution of core co-expressed genes')
#plt.xlabel('Percent of core co-expressed genes')
#plt.ylabel('Density')


# In[ ]:


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
