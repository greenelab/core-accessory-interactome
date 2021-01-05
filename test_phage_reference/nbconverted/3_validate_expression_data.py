
# coding: utf-8

# # Validate expression data

# In[1]:


import seaborn as sns
import pandas as pd
from core_acc_modules import paths


# In[2]:


# Load expression data

expression_pao1_df = pd.read_csv(paths.PAO1_GE, sep='\t', index_col=0, header=0)
expression_pa14_df = pd.read_csv(paths.PA14_GE, sep='\t', index_col=0, header=0)
expression_phage_df = pd.read_csv(paths.PHAGE_GE, sep='\t', index_col=0, header=0)


# In[3]:


expression_phage_df


# In[4]:


expression_pao1_df.head()


# In[5]:


expression_pa14_df.head()


# ## Negative case
# 
# Check that *E. Coli* sample (SRR13234437) poorly aligns with all reference genomes

# In[6]:


# Check there is no expression in E. Coli sample using phage reference
print("Total number of genes:", len(expression_phage_df.loc["SRR13234437"]))
print("Percent of 0 expressed genes:",
      (expression_phage_df.loc["SRR13234437"] == 0).sum()/len(expression_phage_df.loc["SRR13234437"]))
print(expression_phage_df.loc["SRR13234437"].max())
sns.distplot(expression_phage_df.loc["SRR13234437"], kde=False)


# In[7]:


# Check there is no expression in E. Coli sample using PAO1 reference
print("Total number of genes:", len(expression_pao1_df.loc["SRR13234437"]))
print("Percent of 0 expressed genes:", 
      (expression_pao1_df.loc["SRR13234437"] == 0).sum()/len(expression_pao1_df.loc["SRR13234437"]))
print(expression_pao1_df.loc["SRR13234437"].max())
sns.distplot(expression_pao1_df.loc["SRR13234437"], kde=False)


# In[ ]:


# Check there is no expression in E. Coli sample using PA14 reference
print("Total number of genes:", len(expression_pa14_df.loc["SRR13234437"]))
print("Percent of 0 expressed genes:", 
      (expression_pa14_df.loc["SRR13234437"] == 0).sum()/len(expression_pa14_df.loc["SRR13234437"]))
print(expression_pa14_df.loc["SRR13234437"].max())
sns.distplot(expression_pa14_df.loc["SRR13234437"], kde=False)


# ## Positive case: PAO1
# 
# * Check that PAO1 sample (SRR13160334) aligns well to PAO1 and phage reference. 
# * Check that PAO1 sample has expression at LUZ19 (NC_010326.1 in phage genome) and Pf1 genes (PA0717-PA0726 PAO1 genome)
# * PA0717-PA0726 = PGD104183, PGD104185, PGD104187, PGD104189, PGD104191, PGD104193, PGD104195, PGD104197, PGD104199, PGD104201

# In[10]:


# Check there is expression in PAO1 sample using PAO1 reference
print("Total number of genes:", len(expression_pao1_df.loc["SRR13160334"]))
print("Percent of 0 expressed genes:", 
      (expression_pao1_df.loc["SRR13160334"] > 0).sum()/len(expression_pao1_df.loc["SRR13160334"]))
print(expression_pao1_df.loc["SRR13160334"].max())
sns.distplot(expression_pao1_df.loc["SRR13160334"])


# In[13]:


"""# Check there is expression in PAO1 sample using PAO1 reference genes: PA0717-PA0726
pao1_pf1s = [
    "PGD104183",
    "PGD104185",
    "PGD104187",
    "PGD104189",
    "PGD104191",
    "PGD104193", 
    "PGD104195", 
    "PGD104197",
    "PGD104199", 
    "PGD104201"
    ]
expression_pao1_df.loc["SRR13160334", pao1_pf1s]"""


# In[11]:


"""# Check there is expression in PAO1 sample using phage reference
print("Total number of genes:", len(expression_phage_df.loc["SRR13160334"]))
print("Percent of 0 expressed genes:", 
      (expression_phage_df.loc["SRR13160334"] > 0).sum()/len(expression_phage_df.loc["SRR13160334"]))
print(expression_phage_df.loc["SRR13160334"].max())
sns.distplot(expression_phage_df.loc["SRR13160334"], kde=True)"""


# In[12]:


# Check there is expression in PAO1 sample using phage reference gene: NC_010326.1
expression_phage_df.loc["SRR13160334", "NC_010326.1"]


# ## Positive case: PA14
# 
# * Check that PA14 sample (ERR3642743) aligns well to PA14 and phage references. 
# * Check that PA14 sample has expression at Pf5 (PA14_49010, PA14_49020 in PA14 genome) and Pf1 genes (AY324828.1, NC_001331.1, MG250485.1 in phage genome)
# * PA14_49010, PA14_49020 = PGD1658748, PGD1658750

# In[14]:


# Check there is expression in PA14 sample using PA14 reference
print("Total number of genes:", len(expression_pa14_df.loc["ERR3642743"]))
print("Percent of 0 expressed genes:", 
      (expression_pa14_df.loc["ERR3642743"] > 0).sum()/len(expression_pa14_df.loc["ERR3642743"]))
print(expression_pa14_df.loc["ERR3642743"].max())
sns.distplot(expression_pa14_df.loc["ERR3642743"], kde=False)


# In[15]:


# Check there is expression in PA14 sample using PA14 reference genes: PA14_49010, PA14_49020
expression_pa14_df.loc["ERR3642743", ["PGD1658748", "PGD1658750"]]


# In[ ]:


# Check there is expression in PA14 sample using phage reference
#sns.distplot(expression_phage_df.loc["ERR3642743"], kde=False)


# In[ ]:


# Check there is expression in PA14 sample using phage reference genes: 
#expression_phage_df.loc["ERR3642743", ["AY324828.1", "NC_001331.1", "MG250485.1"]]


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
