#!/usr/bin/env python
# coding: utf-8

# # Visualize expression data

# In[1]:


import seaborn as sns
import pandas as pd
from pathlib import Path
from core_acc_modules import utils, paths_phage


# In[2]:


# Load expression data
expression_phage_df = pd.read_csv(paths_phage.PHAGE_GE, sep='\t', index_col=0, header=0)


# In[3]:


# Plot quantification for negative samples
sns.distplot(expression_phage_df.mean(), kde=False)


# In[ ]:


# Plot quantification for positive sample

