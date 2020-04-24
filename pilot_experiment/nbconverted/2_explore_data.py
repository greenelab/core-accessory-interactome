
# coding: utf-8

# # Explore data
# This notebook performs a first pass exploration of the data. In particular, this notebook examines the types of interactions that exist between genes

# In[2]:


import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[3]:


# Input
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

real_expression_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "input",
    "selected_normalized_data.tsv")

shuffled_expression_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "input",
    "shuffled_selected_normalized_data.tsv")

gene_annot_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "annotations",
    "selected_gene_annotations.txt")


# In[11]:


# Read data
real_expression = pd.read_table(
    real_expression_file,
    header=0,
    sep='\t',
    index_col=0)

shuffled_expression = pd.read_table(
    shuffled_expression_file,
    header=0,
    sep='\t',
    index_col=0)

gene_annot = pd.read_table(
    gene_annot_file,
    header=0,
    sep='\t',
    index_col=0)

real_expression.head()


# In[4]:


shuffled_expression.head()


# In[5]:


gene_annot.head()


# In[6]:


# Group genes by core and accessory annotation
core_gene_ids = list(gene_annot[gene_annot['annotation'] == 'core'].index)
acc_gene_ids = list(gene_annot[gene_annot['annotation'] == 'accessory'].index)


# # Examine interactions between genes
# 
# First we can examine how coordinated groups of genes are: core-core, core-accessory, accessory-accessory 

# ## Get groups of expression data

# In[7]:


# Get core-core gene expression
real_core_expression = real_expression[core_gene_ids]


# In[8]:


# Get accessory-accessory gene expression
real_acc_expression = real_expression[acc_gene_ids]


# ## Calculate correlation between gene profiles

# In[9]:


# Get correlation of core-core genes
real_core_corr = real_core_expression.corr(method='pearson')


# In[10]:


# Get correlation of accessory-accessory genes
real_acc_corr = real_acc_expression.corr(method='pearson')


# In[11]:


# Get correlation of all genes
real_all_corr = real_expression.corr(method='pearson')


# In[12]:


# Get correlation of core-accessory genes
real_core_acc_corr = real_all_corr.loc[core_gene_ids, acc_gene_ids]


# In[13]:


# Get correlation of control dataset
shuffled_all_corr = shuffled_expression.corr(method='pearson')


# ## Plot distribution of correlation scores

# In[15]:


# Flatten and get only upper triangle values from correlation matrix
# Core gene correlations
real_core_corr_score = real_core_corr.values[np.triu_indices(n=len(real_core_corr), k=1)]
core_labels = np.repeat(['core'],len(real_core_corr_score))

real_core_corr_df = pd.DataFrame(data={'corr_score': real_core_corr_score,
                                      'group': core_labels})

print(real_core_corr_df.shape)
print('mean ', np.mean(real_core_corr_score))
print('median ', np.median(real_core_corr_score))
print('var ', np.var(real_core_corr_score))
real_core_corr_df.head()


# In[16]:


# Flatten and get only upper triangle values from correlation matrix
# Accessory gene correlations
real_acc_corr_score = real_acc_corr.values[np.triu_indices(n=len(real_acc_corr), k=1)]
acc_labels = np.repeat(['accessory'],len(real_acc_corr_score))

real_acc_corr_df = pd.DataFrame(data={'corr_score': real_acc_corr_score,
                                     'group': acc_labels})

print(real_acc_corr_df.shape)
print('mean ', np.mean(real_acc_corr_score))
print('median ', np.median(real_acc_corr_score))
print('var ', np.var(real_acc_corr_score))
real_acc_corr_df.head()


# In[17]:


# Flatten and get only upper triangle values from correlation matrix
# Core-accessory gene correlations
real_core_acc_corr_score = real_core_acc_corr.values.flatten().tolist()
core_acc_labels = np.repeat(['core-acc'],len(real_core_acc_corr_score))

real_core_acc_corr_df = pd.DataFrame(data={'corr_score': list(real_core_acc_corr_score),
                                           'group': core_acc_labels})

print(real_core_acc_corr_df.shape)
print('mean ', np.mean(real_core_acc_corr_score))
print('median ', np.median(real_core_acc_corr_score))
print('var ', np.var(real_core_acc_corr_score))
real_core_acc_corr_df.head()


# In[18]:


# Flatten and get only upper triangle values from correlation matrix
# All gene correlations
real_all_corr_score = real_all_corr.values[np.triu_indices(n=len(real_all_corr), k=1)]
all_labels = np.repeat(['all'],len(real_all_corr_score))

real_all_corr_df = pd.DataFrame(data={'corr_score': real_all_corr_score,
                                     'group': all_labels})

print(real_all_corr_df.shape)
print('mean ', np.mean(real_all_corr_score))
print('median ', np.median(real_all_corr_score))
print('var ', np.var(real_all_corr_score))
real_all_corr_df.head()


# In[19]:


# Flatten and get only upper triangle values from correlation matrix
# Shuffled gene correlations
shuffled_all_corr_score = shuffled_all_corr.values[np.triu_indices(n=len(shuffled_all_corr), k=1)]
shuffled_all_labels = np.repeat(['shuffled'],len(real_all_corr_score))

shuffled_all_corr_df = pd.DataFrame(data={'corr_score': shuffled_all_corr_score,
                                          'group': shuffled_all_labels})

print(shuffled_all_corr_df.shape)
print('mean ', np.mean(shuffled_all_corr_score))
print('median ', np.median(shuffled_all_corr_score))
print('var ', np.var(shuffled_all_corr_score))
shuffled_all_corr_df.head()


# In[20]:


# Create df
corr_scores_df = pd.concat([real_core_corr_df,
                            real_acc_corr_df,
                            real_core_acc_corr_df,
                            real_all_corr_df,
                            shuffled_all_corr_df],
                           ignore_index=True)

print(corr_scores_df.shape)
corr_scores_df.head()


# In[21]:


# Plot all correlation scores
sns.boxplot(data=corr_scores_df,
           x='group',
           y='corr_score',
           palette='Set3').set_title('Distribution of correlation scores per group')


# In[48]:


# Distribution plot for core genes
sns.distplot(real_core_corr_score, label='core', color='red')
sns.distplot(real_acc_corr_score, label='accessory', color='blue')
sns.distplot(real_core_acc_corr_score, label='core-accessory', color='purple')
sns.distplot(shuffled_all_corr_score, label='shuffled', color='grey')

plt.legend(prop={'size': 12})
plt.title('Density of correlation scores per group')
plt.ylabel('Density')


# **Some sanity checks**
# * Shuffled dataset has very little correlation between genes, as expected since we have disrupted the inherent relationships between genes through our permutation process
# * Since core genes comprise 97% of the total genes, the mean correlation for all genes is the same as the core gene set
# 
# **Other observations**
# * The mean correlation of accessory genes only is 0.2295, which is higher compared to the mean correlation of core genes, 0.01037.
# * The mean correlation of accessory-core genes is -0.036, which is very slightly lower compared to core-only and accessory only groups.
# 
# * Looking at the density plot for the accessory-accessory gene correlation scores, the scores are shifted to the right.
# * Looking at the density plot for accessory-core gene correlation scores, the scores are shifted very slightly to the left possibly indicating an antagonistic relationship. Need to read more about what is currently known about core and accessory genes and about negative regulation patterns.
# * Looking at the density plot for the core-core gene correlation scores, there is a normal distribution of correlation scores which would imply that a minor proportion of genes are highly positively correlated (i.e. genes in operons) or highly negatively correlated (example?), while the majority of genes are not strongly correlated. Need to read more to determine if this makes sense

# ## Binarize the correlation matrix to get interactions

# In[23]:


# Binarize correlation score to get possible interactions
threshold = 0.9
real_core_edges = real_core_corr>threshold
real_acc_edges = real_acc_corr>threshold
real_core_acc_edges = real_core_acc_corr>threshold
real_all_edges = real_all_corr>threshold
shuffled_all_edges = shuffled_all_corr>threshold


# In[25]:


real_core_edges_score = real_core_edges.values[np.triu_indices(n=len(real_core_edges), k=1)]
real_acc_edges_score = real_acc_edges.values[np.triu_indices(n=len(real_acc_edges), k=1)]
real_core_acc_edges_score = real_core_acc_edges.values.flatten().tolist()
real_all_edges_score = real_all_edges.values[np.triu_indices(n=len(real_all_edges), k=1)]
shuffled_all_edges_score = shuffled_all_edges.values[np.triu_indices(n=len(shuffled_all_edges), k=1)]


# In[29]:


num_edges_df = pd.DataFrame(data={'group': ['core', 'accessory', 'core-accessory', 'all', 'shuffled'],
                                  'num_edges': [sum(real_core_edges_score),
                                                sum(real_acc_edges_score),
                                                sum(real_core_acc_edges_score),
                                                sum(real_all_edges_score),
                                                sum(shuffled_all_edges_score)],
                                  'percent_edges': [sum(real_core_edges_score)/real_core_corr_df.shape[0],
                                                    sum(real_acc_edges_score)/real_acc_corr_df.shape[0],
                                                    sum(real_core_acc_edges_score)/((real_core_acc_corr.shape[0]
                                                                                    *real_core_acc_corr.shape[1])),
                                                    sum(real_all_edges_score)/real_all_corr_df.shape[0],
                                                    sum(shuffled_all_edges_score)/shuffled_all_corr_df.shape[0]]
                                 })

num_edges_df.head()

