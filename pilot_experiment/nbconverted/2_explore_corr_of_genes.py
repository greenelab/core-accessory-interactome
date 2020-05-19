
# coding: utf-8

# # Explore correlation of genes
# This notebook performs a first pass exploration of the data. In particular, this notebook examines the types of interactions that exist between genes and how coordinated groups of genes are: core-core, core-accessory, accessory-accessory 

# In[1]:


import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle

np.random.seed(123)


# In[2]:


# User - which experiments to use
which_experiments = "PA14"


# In[3]:


# Input
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

real_expression_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "input",
    "selected_"+which_experiments+"_normalized_data.tsv")

shuffled_expression_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "input",
    "shuffled_"+which_experiments+"_selected_normalized_data.tsv")

gene_annot_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "annotations",
    "selected_gene_annotations.txt")


# In[4]:


# Output directory to store gene ids and correlation outputs
base_intermediate_dir = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "tmp")


# In[5]:


# Read data
real_expression = pd.read_csv(
    real_expression_file,
    header=0,
    sep='\t',
    index_col=0)

shuffled_expression = pd.read_csv(
    shuffled_expression_file,
    header=0,
    sep='\t',
    index_col=0)

gene_annot = pd.read_csv(
    gene_annot_file,
    header=0,
    sep='\t',
    index_col=0)

real_expression.head()


# In[6]:


shuffled_expression.head()


# In[7]:


gene_annot.head()


# In[8]:


# Group genes by core and accessory annotation
core_gene_ids = list(gene_annot[gene_annot['annotation'] == 'core'].index)
acc_gene_ids = list(gene_annot[gene_annot['annotation'] == 'accessory'].index)


# In[9]:


# Pickle gene ids
pickle.dump(core_gene_ids, open(
    os.path.join(
        base_intermediate_dir,
        "core_gene_ids.pickle"),
    "wb"))
pickle.dump(acc_gene_ids, open(
    os.path.join(
        base_intermediate_dir,
        "acc_gene_ids.pickle"),
    "wb"))


# ## Get groups of expression data

# In[10]:


# Get core-core gene expression
real_core_expression = real_expression[core_gene_ids]


# In[11]:


# Get accessory-accessory gene expression
real_acc_expression = real_expression[acc_gene_ids]


# ## Calculate correlation between gene profiles

# In[12]:


# Get correlation of core-core genes
real_core_corr = real_core_expression.corr(method='pearson')


# In[13]:


# Get correlation of accessory-accessory genes
real_acc_corr = real_acc_expression.corr(method='pearson')


# In[14]:


# Get correlation of all genes
real_all_corr = real_expression.corr(method='pearson')

# Save 
pickle.dump(real_all_corr, open(
    os.path.join(
        base_intermediate_dir,
        "real_all_corr.pickle"),
    "wb"))


# In[15]:


# Get correlation of core-accessory genes
real_core_acc_corr = real_all_corr.loc[core_gene_ids, acc_gene_ids]


# In[16]:


# Get correlation of control dataset
shuffled_all_corr = shuffled_expression.corr(method='pearson')

# Save
pickle.dump(shuffled_all_corr, open(
    os.path.join(
        base_intermediate_dir,
        "shuffled_all_corr.pickle"),
    "wb"))


# ## Plot distribution of correlation scores

# In[17]:


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


# In[18]:


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


# In[19]:


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


# In[20]:


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


# In[21]:


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


# In[22]:


# Create df
corr_scores_df = pd.concat([real_core_corr_df,
                            real_acc_corr_df,
                            real_core_acc_corr_df,
                            real_all_corr_df,
                            shuffled_all_corr_df],
                           ignore_index=True)

print(corr_scores_df.shape)
corr_scores_df.head()


# In[23]:


sns.set()


# In[24]:


# Plot all correlation scores
sns.boxplot(data=corr_scores_df,
           x='group',
           y='corr_score',
           palette='Set3').set_title('Distribution of correlation scores per group')


# In[25]:


# Distribution plot for core genes
sns.distplot(real_core_corr_score, 
             label='core', 
             color='red',
             hist = False, 
             kde = True,
             kde_kws = {'shade': True}
            )

sns.distplot(real_acc_corr_score,
             label='accessory',
             color='blue',
             hist = False, 
             kde = True,
             kde_kws = {'shade': True}
            )
sns.distplot(real_core_acc_corr_score, 
             label='core-accessory', 
             color='purple',
             hist = False, 
             kde = True,
             kde_kws = {'shade': True}
            )
sns.distplot(shuffled_all_corr_score,
             label='shuffled', 
             color='grey',
             hist = False, 
             kde = True,
             kde_kws = {'shade': True}
            )

plt.legend(prop={'size': 12})
plt.title('Probability density of correlation scores per group')
plt.ylabel('Probability Density')


# In[26]:


# Get bins using all data
hist, bins_corr = np.histogram(np.concatenate([real_core_corr_score,
                                              real_acc_corr_score,
                                              real_core_acc_corr_score,
                                              shuffled_all_corr_score]))

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(8,8))

# Distribution plot for core genes
sns.distplot(real_core_corr_score, 
             label='core', 
             color='red',
             bins=bins_corr,
             kde=False,
             ax=axes[0,0]
            )

sns.distplot(real_acc_corr_score,
             label='accessory',
             color='blue',
             bins=bins_corr,
             kde=False,
             ax=axes[0,1]
            )
sns.distplot(real_core_acc_corr_score, 
             label='core-accessory', 
             color='purple',
             bins=bins_corr,
             kde=False,
             ax=axes[1,0]
            )
sns.distplot(shuffled_all_corr_score,
             label='shuffled', 
             color='grey',
             bins=bins_corr,
             kde=False,
             ax=axes[1,1]
            )

plt.suptitle('Histogram of correlation scores per group')
axes[0,0].set_title('Core-Core')
axes[0,1].set_title('Accessory-Accessory')
axes[1,0].set_title('Core-Accessory')
axes[1,1].set_title('Shuffled')
fig.text(0.5, 0.01, 'Correlation between genes', ha='center', fontsize=12)
fig.text(0.01, 0.5, 'Count', ha='center', rotation=90, fontsize=12)
plt.tight_layout(pad=0.4, 
                 w_pad=0.5,
                 h_pad=1.0,
                 rect=[0, 0.03, 1, 0.95])


# **Note about visualizations:**
# * Based on the density plot, we observed a shift in the accessory-accessory gene correlation scores. This density plot represents the probability of a random variable falling within a particular range of values (P(0 <= X <= 0.5)). This probability is given by the integral of this variable's PDF over that range -- that is it is given by the area under the density function. So there is a higher likelihood of correlation scores > 0.5 for accessory-accessory genes compared to core-core genes. But the *exact probability* of accessory-accessory genes having a correlation score > 0.5 is not known from just visually inspecting the density plot. To get the exact probability we would need to calculate the integral from 0.5 onward.
# * While the shift if very clear to see in the density plot, the meaning of the y-axis is not as intuitive, so we also plot histograms for each group of genes. The histogram also shows a shift in the accessory-accessory gene correlation scores and also shows that the number of accessory-accessory interactions is orders of magnitude lower compared to core-core interactions.
# 
# https://towardsdatascience.com/histograms-and-density-plots-in-python-f6bda88f5ac0
# 
# **Some sanity checks:**
# * Shuffled dataset has very little correlation between genes, as expected since we have disrupted the inherent relationships between genes through our permutation process
# * Since core genes comprise 97% of the total genes, the mean correlation for all genes is the same as the core gene set
# 
# **Overall observations:**
# * Looking at the density plot for the accessory-accessory gene correlation scores, the scores are shifted to the right.
# * The shift in acc-acc genes having higher correlation scores is more prominent in mixed sample datasets (`which_experiment == 'All'`). This shift is especially prominent in the PA14-only dataset (`which_experiment == 'PA14'`), where the accessory genes are absent.

# ## Examine expression of genes per group
# **Question**
# Is the reason for this shift because the accessory genes are absent = expression values for those accessory genes are all very low with little variance?

# In[27]:


# Get mean, max, min expression per core gene
mean_real_core_expression = real_core_expression.mean()
max_real_core_expression = real_core_expression.max()
min_real_core_expression = real_core_expression.min()


# In[28]:


# Get mean, max, min expression per accessory gene
mean_real_acc_expression = real_acc_expression.mean()
max_real_acc_expression = real_acc_expression.max()
min_real_acc_expression = real_acc_expression.min()


# In[29]:


# Distribution plot for core genes
sns.distplot(mean_real_core_expression, 
             label='core', 
             color='red',
             hist = False, 
             kde = True,
             kde_kws = {'shade': True}
            )

sns.distplot(mean_real_acc_expression,
             label='accessory',
             color='blue',
             hist = False, 
             kde = True,
             kde_kws = {'shade': True}
            )

plt.legend(prop={'size': 12})
plt.title('Probability density of mean gene expression')
plt.ylabel('Probability Density')
plt.xlabel('Mean gene expression')


# In[30]:


# Get bins using all data
hist, bins_expression = np.histogram(np.concatenate([mean_real_core_expression,
                                              mean_real_acc_expression]))


# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(6,6))

# Distribution plot for core genes
sns.distplot(mean_real_core_expression, 
             label='core', 
             color='red',
             bins=bins_expression,
             kde=False,
             ax=axes[0]
            )

sns.distplot(mean_real_acc_expression,
             label='accessory',
             color='blue',
             bins=bins_expression,
             kde=False,
             ax=axes[1]
            )

plt.suptitle('Histogram of gene expression per group')
axes[0].set_title('Core-Core')
axes[1].set_title('Accessory-Accessory')
fig.text(0.5, 0.01, 'Gene expression', ha='center', fontsize=12)
fig.text(0.01, 0.5, 'Count', ha='center', rotation=90, fontsize=12)
plt.tight_layout(pad=0.4, 
                 w_pad=0.5,
                 h_pad=1.0,
                 rect=[0, 0.03, 1, 0.95])

