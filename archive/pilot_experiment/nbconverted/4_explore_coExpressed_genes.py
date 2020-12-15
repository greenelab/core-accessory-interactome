
# coding: utf-8

# # Explore highly co-expressed genes
# In the previous [notebook](2_explore_corr_of_genes.ipynb) we observed that using 39 samples with 201 PAO1-specific genes, that the correlation of accessory-accessory genes is higher compared to the correlation of core-core and core-accessory genes.
# 
# Based on this finding, we want to know: *What can explain this difference in correlation distribution?*
# 
# This notebook performs a follow-up analysis. In particular this notebook performs a deeper examination of the correlation structure per group (core-core, core-accessory, accessory-accessory) by looking at the trends of the nearest neighbors (i.e. highly correlated genes) of each gene.

# In[74]:


import pandas as pd
import os
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from functions import calculations

np.random.seed(123)


# In[75]:


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
operon_file = "https://github.com/greenelab/adage/blob/master/Genome_organization/operon_3.txt"


# In[76]:


# Read in gene ids
core_gene_ids = pickle.load(open(core_gene_ids_file, "rb"))
acc_gene_ids = pickle.load(open(acc_gene_ids_file, "rb"))

# Get number of core and accessory genes
num_core_genes = len(core_gene_ids)
num_acc_genes = len(acc_gene_ids)
num_all_genes = num_core_genes + num_acc_genes


# # Extract statistics about co-expression from correlation matrix

# In[77]:


# Define threshold for highly co-expressed genes
coexpression_threshold = 0.9


# ### Co-expression statsitics using real data

# In[78]:


# Get co-expression patterns using real expression data
real_core_df, real_acc_df = calculations.get_coexpression_stats(real_all_corr_file,
                                                                operon_file,
                                                                core_gene_ids_file,
                                                                acc_gene_ids_file,
                                                                coexpression_threshold)


# In[79]:


real_core_df.head()


# In[80]:


real_acc_df.head()


# ### Co-expression statistics using shuffled data

# In[81]:


# Get co-expression patterns using shuffled expression data (control)
shuffled_core_df, shuffled_acc_df = calculations.get_coexpression_stats(shuffled_all_corr_file,
                                                                operon_file,
                                                                core_gene_ids_file,
                                                                acc_gene_ids_file,
                                                                coexpression_threshold)

shuffled_core_df.head()


# In[82]:


shuffled_acc_df.head()


# # Plot trends in co-expression data

# ## 1. Number of co-expressed genes

# In[83]:


sns.set()


# ### Number of co-expressed genes

# In[84]:


# Get bins using all data
hist, bins_num_coexpressed_real = np.histogram(np.concatenate([real_core_df['num_coexpressed_genes'].values,
                                                               real_acc_df['num_coexpressed_genes'].values]))

# Distribution of number of co-expressed genes in real data
# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

sns.distplot(real_core_df['num_coexpressed_genes'].tolist(), 
             label='core',
             color='red',
             kde=False,
             bins=bins_num_coexpressed_real,
             ax=axes[0])

sns.distplot(real_acc_df['num_coexpressed_genes'].tolist(), 
             label='accessory',
             color='blue',
             kde=False,
             bins=bins_num_coexpressed_real,
             ax=axes[1])

plt.suptitle('Number of co-expressed genes (real data, threshold={})'.format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Number of co-expressed genes', ha='center')
axes[0].set_ylabel('Counts')


# ### Number of nonzero co-expressed genes

# In[85]:


## Remove genes with 0 co-expressed genes
# Get bins using all data
hist, bins_num_coexpressed_real_nonzero = np.histogram(np.concatenate(
    [real_core_df[real_core_df['num_coexpressed_genes']>0]['num_coexpressed_genes'].values,
     real_acc_df[real_acc_df['num_coexpressed_genes']>0]['num_coexpressed_genes'].values]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution of number of co-expressed genes in real data
sns.distplot(real_core_df[real_core_df['num_coexpressed_genes']>0]['num_coexpressed_genes'].tolist(), 
             label='core',
             color='red',
             kde=False,
             bins=bins_num_coexpressed_real_nonzero,
             ax=axes[0])

sns.distplot(real_acc_df[real_acc_df['num_coexpressed_genes']>0]['num_coexpressed_genes'].tolist(), 
             label='accessory',
             color='blue',
             kde=False,
             bins=bins_num_coexpressed_real_nonzero,
             ax=axes[1])

plt.suptitle('Number of nonzero co-expressed genes (real data, threshold={})'.format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Number of co-expressed genes', ha='center')
axes[0].set_ylabel('Counts')


# In[86]:


# Get bins using all data
hist, bins_num_coexpressed_shuffled = np.histogram(np.concatenate([shuffled_core_df['num_coexpressed_genes'].values,
                                                                   shuffled_acc_df['num_coexpressed_genes'].values]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution of number of co-expressed genes in shuffled data
sns.distplot(shuffled_core_df['num_coexpressed_genes'].tolist(), 
             label='core',
             color='red',
             kde=False,
             bins=bins_num_coexpressed_shuffled,
             ax=axes[0])

sns.distplot(shuffled_acc_df['num_coexpressed_genes'].tolist(), 
             label='accessory', 
             color='blue',
             kde=False,
             bins=bins_num_coexpressed_shuffled,
             ax=axes[1]
            )

plt.suptitle('Number of co-expressed genes (shuffled data, threshold={})'.format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Number of co-expressed genes', ha='center')
axes[0].set_ylabel('Counts')


# In[87]:


# Print statistics about co-expressed genes
print('Using a threshold of {} to define co-expression (real data): \n'.
     format(coexpression_threshold))
print('- For a given CORE gene, there is a median of {} co-expressed  genes'.
      format(np.median(real_core_df['num_coexpressed_genes'])))
print('- For a given ACCESSORY gene, there is a median of {} co-expressed  genes \n'.
      format(np.median(real_acc_df['num_coexpressed_genes'])))

# For shuffled data
print('Using a threshold of {} to define co-expression (shuffled data): \n'.
     format(coexpression_threshold))
print('- For a given CORE gene, there is a median of {} co-expressed  genes'.
      format(np.median(shuffled_core_df['num_coexpressed_genes'])))
print('- For a given ACCESSORY gene, there is a median of {} co-expressed  genes'.
      format(np.median(shuffled_acc_df['num_coexpressed_genes'])))


# **Overall:**
# * Many core and accessory genes are not co-expressed with other genes or very few genes, as expected
# * As increase threshold (more stringent), there are fewer co-expressed genes, as expected
# * (control, not shown) All genes (threshold=0.75,0.9) are connected to just 1 gene, as expected, since we have destroyed relationships between genes when we shuffled
# * Both core and accessory genes are have the same degree of connectivity
# 
# **Observation using All experiments:**
# * At a threshold of 0.5, core genes have a median of 716 genes they are connect to and accessory genes have a median of 601 genes they are connected to
# * At a threshold of 0.75 core genes have a median of 53 genes they are connect to and accessory genes have a median of 74 genes they are connected to
# * At a threshold of 0.9, core and accessory genes have a median of 1 gene that they are connected to 
# 
# **Observation using only PAO1 experiments:**
# * At a threshold of 0.5, core genes have a median of 1331 genes they are connect to and accessory genes have a median of 1341 genes they are connected to. So there are more connections using only PAO1 genes, which we would expect since the accessory genes are PAO1-specific.
# * At a threshold of 0.75 core genes have a median of 554 genes they are connect to and accessory genes have a median of 576 genes they are connected to
# * At a threshold of 0.9, core genes have a median of 21 gene that they are connected to and accessory genes have a median of 57 genes they are connected to 
# 
# **Observation using only PA14 experiments:**
# * At a threshold of 0.5, core genes have a median of 1133 genes they are connect to and accessory genes have a median of 1065 genes they are connected to. So there are more connections using only PAO1 genes, which we would expect since the accessory genes are PAO1-specific.
# * At a threshold of 0.75 core genes have a median of 153 genes they are connect to and accessory genes have a median of 132 genes they are connected to
# * At a threshold of 0.9, core genes have a median of 2 gene that they are connected to and accessory genes have a median of 3 genes they are connected to 

# ### Compare co-expressed vs not co-expressed genes
# What is the difference between genes that are co-expressed vs those that are not?

# In[15]:


# Get genes that are co-expressed with other genes
coexpressed_core_genes = list(real_core_df[real_core_df['num_coexpressed_genes']>0].index)
coexpressed_acc_genes = list(real_acc_df[real_acc_df['num_coexpressed_genes']>0].index)


# In[16]:


len(coexpressed_core_genes)


# In[17]:


len(coexpressed_acc_genes)


# At a strict threshold of 0.9, all genes are co-expressed with **at least** one other gene. So there are no independent genes.

# ## 2. Percent of co-expressed genes that are NOT in the same operon

# In[15]:


# Calculate the percent of co-expressed genes that are non co-operonic (real data)
real_percent_non_cooperonic_coexpressed_core_genes = (
    real_core_df['num_non_cooperonic_coexpressed_genes']/real_core_df['num_coexpressed_genes'])

real_percent_non_cooperonic_coexpressed_acc_genes = (
    real_acc_df['num_non_cooperonic_coexpressed_genes']/real_acc_df['num_coexpressed_genes'])

# There are NaNs in cases where there are 0 co-expressed genes and therefore 0 non-cooperonic genes
real_num_core_na = real_percent_non_cooperonic_coexpressed_core_genes.isna().sum()
real_num_acc_na = real_percent_non_cooperonic_coexpressed_acc_genes.isna().sum()

# Since we are concerned with "of those co-expressed genes how many are in NOT in the same operon", we will remove these
real_percent_non_cooperonic_coexpressed_core_genes_noNa = real_percent_non_cooperonic_coexpressed_core_genes.dropna(inplace=False)
real_percent_non_cooperonic_coexpressed_acc_genes_noNa = real_percent_non_cooperonic_coexpressed_acc_genes.dropna(inplace=False)


# In[16]:


## TEST: What does distribution look like before removing NaNs?

# Get bins using all data
hist, bins_num_percent_non_cooperonic_real = np.histogram(
    np.concatenate(
        [real_percent_non_cooperonic_coexpressed_core_genes.fillna(0),
         real_percent_non_cooperonic_coexpressed_acc_genes.fillna(0)]
    )
)

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution of percent of co-expressed genes that are NOT co-operonic in real data
sns.distplot(real_percent_non_cooperonic_coexpressed_core_genes.fillna(0),
             label='core',
             color='red',
             kde=False,
             bins=bins_num_percent_non_cooperonic_real,
             ax=axes[0])

sns.distplot(real_percent_non_cooperonic_coexpressed_acc_genes.fillna(0), 
             label='accessory',
             color='blue',
             kde=False, 
             bins=bins_num_percent_non_cooperonic_real,
             ax=axes[1])

plt.suptitle('TEST'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes', ha='center')
axes[0].set_ylabel('Count')


# In[17]:


# Get bins using all data
hist, bins_num_percent_non_cooperonic_real = np.histogram(
    np.concatenate(
        [real_percent_non_cooperonic_coexpressed_core_genes_noNa,
         real_percent_non_cooperonic_coexpressed_acc_genes_noNa]
    )
)

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution of percent of co-expressed genes that are NOT co-operonic in real data
sns.distplot(real_percent_non_cooperonic_coexpressed_core_genes_noNa,
             label='core',
             color='red',
             kde=False,
             bins=bins_num_percent_non_cooperonic_real,
             ax=axes[0])

sns.distplot(real_percent_non_cooperonic_coexpressed_acc_genes_noNa, 
             label='accessory',
             color='blue', 
             kde=False,
             bins=bins_num_percent_non_cooperonic_real,
             ax=axes[1]
            )

plt.suptitle('Percent of co-expressed genes that are NOT co-operonic (real data, threshold={})'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes', ha='center')
axes[0].set_ylabel('Count')


# In[18]:


# Calculate the percent of co-expressed genes that are non co-operonic (shuffled data)
shuffled_percent_non_cooperonic_coexpressed_core_genes = (
    shuffled_core_df['num_non_cooperonic_coexpressed_genes']/shuffled_core_df['num_coexpressed_genes'])

shuffled_percent_non_cooperonic_coexpressed_acc_genes = (
    shuffled_acc_df['num_non_cooperonic_coexpressed_genes']/shuffled_acc_df['num_coexpressed_genes'])

# There are NaNs in cases where there are 0 co-expressed genes and therefore 0 non-cooperonic genes
shuffled_num_core_na = shuffled_percent_non_cooperonic_coexpressed_core_genes.isna().sum()
shuffled_num_acc_na = shuffled_percent_non_cooperonic_coexpressed_acc_genes.isna().sum()

# Since we are concerned with "of those co-expressed genes how many are in NOT in the same operon", we will remove these
shuffled_percent_non_cooperonic_coexpressed_core_genes_noNa = shuffled_percent_non_cooperonic_coexpressed_core_genes.dropna(
    inplace=False)
shuffled_percent_non_cooperonic_coexpressed_acc_genes_noNa = shuffled_percent_non_cooperonic_coexpressed_acc_genes.dropna(
    inplace=False)


# In[19]:


# Get bins using all data
hist, bins_num_percent_non_cooperonic_shuffled = np.histogram(
    np.concatenate([shuffled_percent_non_cooperonic_coexpressed_core_genes_noNa,
                    shuffled_percent_non_cooperonic_coexpressed_acc_genes_noNa]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution of percent of co-expressed genes that are NOT co-operonic in shuffled data
sns.distplot(shuffled_percent_non_cooperonic_coexpressed_core_genes_noNa,
             label='core',
             color='red',
             kde=False,
             bins=bins_num_percent_non_cooperonic_shuffled,
             ax=axes[0])

sns.distplot(shuffled_percent_non_cooperonic_coexpressed_acc_genes_noNa, 
             label='accessory',
             color='blue', 
             kde=False,
             bins=bins_num_percent_non_cooperonic_shuffled,
             ax=axes[1])

plt.suptitle('Percent of co-expressed genes that are NOT co-operonic (shuffled data, threshold={})'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes', ha='center')
axes[0].set_ylabel('Count')


# In[20]:


# Print statistics about non co-operonic co-expressed genes
print('Using a threshold of {} to define co-expression (real data): \n'.
     format(coexpression_threshold))
print('- We removed {} ({}%) genes that had 0 co-expressed genes'.
      format(real_num_core_na, round((real_num_core_na/num_core_genes)*100,2)))
print('''- Of those remaining genes, for a given CORE gene, 
      {}% (median) of the co-expressed genes are NOT in a shared operon'''.
      format(round(np.median(real_percent_non_cooperonic_coexpressed_core_genes_noNa)*100,2)))

print('- We removed {} ({}%) genes that had 0 co-expressed genes'.
      format(real_num_acc_na, round(real_num_acc_na/num_acc_genes)*100,2))
print('''- Of those remaining genes, for a given ACCESSORY gene, 
     {}% (median) of the co-expressed genes are NOT in a shared operon \n'''.
      format(round(np.median(real_percent_non_cooperonic_coexpressed_acc_genes_noNa)*100,2)))


# For shuffled data
print('Using a threshold of {} to define co-expression (shuffled data): \n'.
     format(coexpression_threshold))
print('- We removed {} ({}%) genes that had 0 co-expressed genes'.
      format(shuffled_num_core_na, round((shuffled_num_core_na/num_core_genes)*100,2)))
print('''- Of those remaining genes, for a given CORE gene, 
      {}% (median) of the co-expressed genes are NOT in a shared operon'''.
      format(round(np.median(shuffled_percent_non_cooperonic_coexpressed_core_genes_noNa)*100,2)))

print('- We removed {} ({}%) genes that had 0 co-expressed genes'.
      format(shuffled_num_acc_na, round((shuffled_num_acc_na/num_acc_genes)*100,2)))
print('''- Of those remaining genes, for a given ACCESSORY gene, 
     {}% (median) of the co-expressed genes are NOT in a shared operon'''.
      format(round(np.median(shuffled_percent_non_cooperonic_coexpressed_acc_genes_noNa)*100,2)))


# **Observations:**
# * For majority of core genes, their co-expressed genes are not in the same operon. For example if core gene A is co-expressed with genes {X,Y,Z} and A is in operon {A,B,C} there are no intersecting genes. 
# * This is the case for most core and accessory genes. There are few co-expression gene sets that overlap with known operons.
# * I would've expected more co-expressed genes to overlap with operons but need to read more about *P. aeruginosa* gene-gene interactions to get more of an intuition.
# * Good to see that as you decrease the threshold there are more non co-operonic gene co-expressed gene sets.
# * (control) There are 0 non-cooperonic co-expressed genes since there are 0 co-expressed genes

# ### Core gene relationships
# For a given core gene, there exists a set of genes that are co-expressed with it. What percent of those co-expressed genes are core?
# Similarly, for a given accessory gene, there exists a set of genes that are co-expressed with it. What percent of those co-expressed genes are core? 

# In[21]:


# We only want to consider those genes with some co-expressed genes
real_core_nan_ids = list(
    real_percent_non_cooperonic_coexpressed_core_genes[real_percent_non_cooperonic_coexpressed_core_genes.isna()].index)
real_acc_nan_ids = list(
    real_percent_non_cooperonic_coexpressed_acc_genes[real_percent_non_cooperonic_coexpressed_acc_genes.isna()].index)

# We also only want to consider only those with some non-cooperonic genes
real_core_zero_ids = list(real_core_df[real_core_df['num_non_cooperonic_coexpressed_genes'] == 0].index)
real_acc_zero_ids = list(real_acc_df[real_acc_df['num_non_cooperonic_coexpressed_genes'] == 0].index)

# Get the union of excluded gene ids
real_exclude_core_ids = set(real_core_zero_ids).union(real_core_nan_ids)
real_exclude_acc_ids = set(real_acc_zero_ids).union(real_acc_nan_ids)

# Since we are concerned with "of those co-expressed gene NOT in the same operon", we will remove the above ids
real_percent_core_with_refcore = real_core_df['percent_non_cooperonic_coexpressed_core'].drop(labels=real_exclude_core_ids)
real_percent_core_with_refacc = real_acc_df['percent_non_cooperonic_coexpressed_core'].drop(labels=real_exclude_acc_ids)


# In[22]:


## Test: Distribution without removing 0 expressed and 0 non co-operonic genes

# Get bins using all data
hist, bins_core_real = np.histogram(np.concatenate([real_core_df['percent_non_cooperonic_coexpressed_core'].values,
                                                    real_acc_df['percent_non_cooperonic_coexpressed_core'].values]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution plot for percent of core co-expressed genes in real data
sns.distplot(real_core_df['percent_non_cooperonic_coexpressed_core'].tolist(),
             label='core', 
             color='red',
             kde=False,
             bins=bins_core_real,
             ax=axes[0])

sns.distplot(real_acc_df['percent_non_cooperonic_coexpressed_core'].tolist(),
             label='accessory',
             color='blue',
             kde=False,
             bins=bins_core_real,
             ax=axes[1])

plt.suptitle('TEST'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes that are core', ha='center')
axes[0].set_ylabel('Count')


# In[23]:


# Get bins using all data
hist, bins_core_real = np.histogram(np.concatenate([real_percent_core_with_refcore.values,
                                                    real_percent_core_with_refacc.values]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution plot for percent of core co-expressed genes in real data
sns.distplot(real_percent_core_with_refcore.tolist(),
             label='core', 
             color='red',
             kde=False,
             bins=bins_core_real,
             ax=axes[0])

sns.distplot(real_percent_core_with_refacc.tolist(),
             label='accessory',
             color='blue',
             kde=False,
             bins=bins_core_real,
             ax=axes[1])

plt.suptitle('Percent of co-expressed non-cooperonic genes that are core (real data, threshold={})'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes that are core', ha='center')
axes[0].set_ylabel('Count')


# In[24]:


# We only want to consider those genes with some co-expressed genes
shuffled_core_nan_ids = list(
    shuffled_percent_non_cooperonic_coexpressed_core_genes[
        shuffled_percent_non_cooperonic_coexpressed_core_genes.isna()].index)
shuffled_acc_nan_ids = list(
    shuffled_percent_non_cooperonic_coexpressed_acc_genes[
        shuffled_percent_non_cooperonic_coexpressed_acc_genes.isna()].index)

# We also only want to consider only those with some non-cooperonic genes
shuffled_core_zero_ids = list(shuffled_core_df[shuffled_core_df['num_non_cooperonic_coexpressed_genes'] == 0].index)
shuffled_acc_zero_ids = list(shuffled_acc_df[shuffled_acc_df['num_non_cooperonic_coexpressed_genes'] == 0].index)

# Get the union of excluded gene ids
shuffled_exclude_core_ids = set(shuffled_core_zero_ids).union(shuffled_core_nan_ids)
shuffled_exclude_acc_ids = set(shuffled_acc_zero_ids).union(shuffled_acc_nan_ids)

# Since we are concerned with "of those co-expressed gene NOT in the same operon", we will remove the above ids
shuffled_percent_core_with_refcore = shuffled_core_df['percent_non_cooperonic_coexpressed_core'].drop(labels=shuffled_exclude_core_ids)
shuffled_percent_core_with_refacc = shuffled_acc_df['percent_non_cooperonic_coexpressed_core'].drop(labels=shuffled_exclude_acc_ids)


# In[25]:


# Get bins using all data
hist, bins_core_shuffled = np.histogram(
    np.concatenate([shuffled_percent_core_with_refcore.values,
                    shuffled_percent_core_with_refacc.values]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution plot for percent of core co-expressed genes in shuffled data
sns.distplot(shuffled_percent_core_with_refcore.tolist(), 
             label='core', 
             color='red',
             kde=False,
             bins=bins_core_shuffled,
             ax=axes[0])

sns.distplot(shuffled_percent_core_with_refacc.tolist(),
             label='accessory', 
             color='blue',
             kde=False,
             bins=bins_core_shuffled,
             ax=axes[1])

plt.suptitle('Percent of co-expressed non-cooperonic genes that are core (shuffled data, threshold={})'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes that are core', ha='center')
axes[0].set_ylabel('Count')


# In[26]:


# Print statistics about non co-operonic co-expressed core genes
print('Using a threshold of {} to define co-expression (real data): \n'.
     format(coexpression_threshold))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(real_exclude_core_ids), round((len(real_exclude_core_ids)/num_core_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given CORE gene,
      {}% (median) of co-expressed non co-operonic genes are core genes'''.
      format(round(np.median(real_percent_core_with_refcore)*100,2)))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(real_exclude_acc_ids), round((len(real_exclude_acc_ids)/num_acc_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given ACCESSORY gene,
      {}% (median) of co-expressed non co-operonic genes are core genes'''.
      format(round(np.median(real_percent_core_with_refacc)*100,2)))

# shuffled data
print('Using a threshold of {} to define co-expression (real data): \n'.
     format(coexpression_threshold))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(shuffled_exclude_core_ids), round((len(shuffled_exclude_core_ids)/num_core_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given CORE gene,
      {}% (median) of co-expressed non co-operonic genes are core genes'''.
      format(round(np.median(shuffled_percent_core_with_refcore)*100,2)))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(shuffled_exclude_acc_ids), round((len(shuffled_exclude_acc_ids)/num_acc_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given ACCESSORY gene,
      {}% (median) of co-expressed non co-operonic genes are core genes'''.
      format(round(np.median(shuffled_percent_core_with_refacc)*100,2)))


# ### Accessory gene relationships
# For a given core gene, there exists a set of genes that are co-expressed with it. What percent of those co-expressed genes are accessory?
# Similarly, for a given accessory gene, there exists a set of genes that are co-expressed with it. What percent of those co-expressed genes are accessory? 

# In[27]:


# Since we are concerned with "of those co-expressed gene NOT in the same operon", we will remove the above ids
real_percent_acc_with_refcore = real_core_df['percent_non_cooperonic_coexpressed_acc'].drop(labels=real_exclude_core_ids)
real_percent_acc_with_refacc = real_acc_df['percent_non_cooperonic_coexpressed_acc'].drop(labels=real_exclude_acc_ids)


# In[28]:


# Get bins using all data
hist, bins_acc_real = np.histogram(np.concatenate([real_percent_acc_with_refcore.values,
                                                   real_percent_acc_with_refacc.values]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution plot for percent of accessory co-expressed genes in real data
sns.distplot(real_percent_acc_with_refcore.tolist(),
             label='core', 
             color='red', 
             kde=False,
             bins=bins_acc_real,
             ax=axes[0])

sns.distplot(real_percent_acc_with_refacc.tolist(),
             label='accessory',
             color='blue', 
             kde=False,
             bins=bins_acc_real,
             ax=axes[1])

plt.suptitle('Percent of co-expressed non-cooperonic genes that are accessory (real data, threshold={})'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes that are accessory', ha='center')
axes[0].set_ylabel('Count')


# In[29]:


# Since we are concerned with "of those co-expressed gene NOT in the same operon", we will remove the above ids
shuffled_percent_acc_with_refcore = shuffled_core_df['percent_non_cooperonic_coexpressed_acc'].drop(labels=shuffled_exclude_core_ids)
shuffled_percent_acc_with_refacc = shuffled_acc_df['percent_non_cooperonic_coexpressed_acc'].drop(labels=shuffled_exclude_acc_ids)


# In[30]:


# Get bins using all data
hist, bins_acc_shuffled = np.histogram(
    np.concatenate([shuffled_percent_acc_with_refcore.values,
                    shuffled_percent_acc_with_refacc.values]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1)

# Distribution plot for percent of accessory co-expressed genes in shuffled data
sns.distplot(shuffled_percent_acc_with_refcore.tolist(), 
             label='core', 
             color='red',
             kde=False,
             bins=bins_acc_shuffled,
             ax=axes[0])
sns.distplot(shuffled_percent_acc_with_refacc.tolist(), 
             label='accessory', 
             color='blue', 
             kde=False,
             bins=bins_acc_shuffled,
             ax=axes[1])

plt.suptitle('Percent of co-expressed non-cooperonic genes that are accessory (shuffled data, threshold={})'.
          format(coexpression_threshold))
axes[0].legend(prop={'size': 12})
axes[1].legend(prop={'size': 12})
fig.text(0.5, 0.01, 'Percent of co-expressed non-cooperonic genes that are accessory', ha='center')
axes[0].set_ylabel('Count')


# In[31]:


# Print statistics about non co-operonic co-expressed accessory genes
print('Using a threshold of {} to define co-expression (real data): \n'.
     format(coexpression_threshold))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(real_exclude_core_ids), round((len(real_exclude_core_ids)/num_core_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given CORE gene,
      {}% (median) of co-expressed non co-operonic genes are accessory genes'''.
      format(round(np.median(real_percent_acc_with_refcore)*100,2)))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(real_exclude_acc_ids), round((len(real_exclude_acc_ids)/num_acc_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given ACCESSORY gene,
      {}% (median) of co-expressed non co-operonic genes are accessory genes'''.
      format(round(np.median(real_percent_acc_with_refacc)*100,2)))

# shuffled data
print('Using a threshold of {} to define co-expression (real data): \n'.
     format(coexpression_threshold))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(shuffled_exclude_core_ids), round((len(shuffled_exclude_core_ids)/num_core_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given CORE gene,
      {}% (median) of co-expressed non co-operonic genes are accessory genes'''.
      format(round(np.median(shuffled_percent_acc_with_refcore)*100,2)))
print('- We removed {} ({}%) genes that had 0 co-expressed genes or 0 non-cooperonic genes'.
      format(len(shuffled_exclude_acc_ids), round((len(shuffled_exclude_acc_ids)/num_acc_genes)*100,2)))
print('''- Of the non-coperonic co-expressed genes, for a given ACCESSORY gene,
      {}% (median) of co-expressed non co-operonic genes are accessory genes'''.
      format(round(np.median(shuffled_percent_acc_with_refacc)*100,2)))


# **Observation:**
# * Core genes tend to be co-expressed with **only** other core genes
# * Many accessory genes tend to be co-expressed with **only** other core genes, 
# * Some accessory genes are co-expressed with mix of accessory genes, core genes
# * (control) There are no relationships between core-core, core-accessory, accessory-accessory

# # Manually examine co-expressed and co-operonic genes

# In[32]:


# Get genes where the number of co-expressed genes > 0 AND number of non-operonic genes is 0 
# Meaning that all the co-expressed genes are in the same operon
real_refcore_cooperonic = real_core_df[
    (real_core_df['num_coexpressed_genes'] > 0) & (real_core_df['num_non_cooperonic_coexpressed_genes'] == 0)]

print(real_refcore_cooperonic.shape)
real_refcore_cooperonic.head(10)


# In[33]:


real_refacc_cooperonic = real_acc_df[
    (real_acc_df['num_coexpressed_genes'] > 0) & (real_acc_df['num_non_cooperonic_coexpressed_genes'] == 0)]

print(real_refacc_cooperonic.shape)
real_refacc_cooperonic.head(10)


# In[34]:


# Manually select core reference gene and co-expressed gene set that is 100% co-operonic
# Read in correlation matrix
real_all_corr = pickle.load(open(real_all_corr_file, "rb"))

real_all_corr[real_all_corr.loc['PA1216']>coexpression_threshold]['PA1216']


# Manually selected PA1216. PA1216 is in operon containing PA1216 - PA1221 genes (http://3.209.27.103/feature/show/?id=105200&view=operons)
# 

# In[35]:


# Get genes where the number of co-expressed genes > 0 AND number of non-operonic genes > 0 
# Meaning that all the co-expressed genes are in the same operon
real_refcore_cooperonic = real_core_df[
    (real_core_df['num_coexpressed_genes'] > 0) & (real_core_df['num_non_cooperonic_coexpressed_genes'] > 0)]

print(real_refcore_cooperonic.shape)
real_refcore_cooperonic.head(10)


# In[36]:


real_refacc_cooperonic = real_acc_df[
    (real_acc_df['num_coexpressed_genes'] > 0) & (real_acc_df['num_non_cooperonic_coexpressed_genes'] > 0)]

print(real_refacc_cooperonic.shape)
real_refacc_cooperonic.head(10)


# In[37]:


# what does it say in literature, should they be in operon
# Manually select core reference gene and co-expressed gene set that is 100% non co-operonic

real_all_corr.loc['PA0001', real_all_corr.loc['PA0001']>coexpression_threshold]


# Manually selected PA0001. PA0001 is in operon containing PA0002-PA0004 genes (http://www.pseudomonas.com/feature/show/?id=134012&view=operons), which are **not** co-expressed with PA0001. Why is that?
