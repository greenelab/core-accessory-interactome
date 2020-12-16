
# coding: utf-8

# # Validate expression data

# In[1]:


import seaborn as sns
import pandas as pd
from core_acc_modules import paths


# In[2]:


# Load expression data
expression_phage_df = pd.read_csv(paths.PHAGE_GE, sep='\t', index_col=0, header=0)
expression_pao1_df = pd.read_csv(paths.PAO1_GE, sep='\t', index_col=0, header=0)


# In[21]:


expression_phage_df


# In[22]:


expression_pao1_df.head()


# ## Negative case

# In[16]:


# Check there is no expression in E. Coli sample using phage reference
sns.distplot(expression_phage_df.loc["SRR12922100"], kde=False)


# In[17]:


# Check there is no expression in E. Coli sample using PAO1 reference
sns.distplot(expression_pao1_df.loc["SRR12922100"], kde=False)


# ## Positive case

# In[14]:


# Check there is expression in PAO1 sample using PAO1 reference(MT118305.1)
sns.distplot(expression_pao1_df.loc["SRR11809598"], kde=False)


# In[8]:


# Check there is expression in Epa41 phage sample using phage genome (Epa41 MT118305.1)
expression_phage_df.loc[["SRR13196071", "SRR13196068"], "MT118305.1"]


# In[10]:


expression_phage_df.loc[["SRR13196071", "SRR13196068"]]


# In[18]:


sns.distplot(expression_pao1_df.loc["SRR13196071"])


# In[19]:


sns.distplot(expression_pao1_df.loc["SRR13196068"])


# In[9]:


# Check there is expression in Epa14 phage sample using phage genome (Epa14 NC_050144.1)
expression_phage_df.loc[["SRR13196071", "SRR13196068"],"NC_050144.1"]


# **Initial Takeaways:**
# 
# As expected:
# * E. Coli samples had 0 expression using phage and PAO1 references
# * PAO1 sample had expression using PAO1 reference
# * Phage samples had 0 expression using PAO1 reference
# 
# Issues:
# * Pseudomonas phage Epa14, Epa41 samples had 0 expression using phage reference → Not sure why
# 
