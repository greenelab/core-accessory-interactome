
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


# In[27]:


EColi_phage_ref_series = expression_phage_df.loc["SRR13234437"]
nonzero_genes_EColi_phage_ref = EColi_phage_ref_series.iloc[EColi_phage_ref_series.nonzero()]
print(nonzero_genes_EColi_phage_ref)


# In[42]:


# Check there is no expression in E. Coli sample using PAO1 reference
print("Total number of genes:", len(expression_pao1_df.loc["SRR13234437"]))
print("Percent of 0 expressed genes:", 
      (expression_pao1_df.loc["SRR13234437"] == 0).sum()/len(expression_pao1_df.loc["SRR13234437"]))

print("Percent of 0 expressed genes:", 
      ((expression_pao1_df.loc["SRR13234437"] > 0) & (expression_pao1_df.loc["SRR13234437"] < 1000))
       .sum()/len(expression_pao1_df.loc["SRR13234437"]))


# In[48]:


EColi_pao1_ref_series = expression_pao1_df.loc["SRR13234437"]
nonzero_genes_EColi_pao1_ref = EColi_pao1_ref_series.iloc[EColi_pao1_ref_series.nonzero()]
sns.distplot(nonzero_genes_EColi_pao1_ref, kde=False)

# Get gene with the max expression
EColi_pao1_ref_series[EColi_pao1_ref_series > 290000]


# In[41]:


# Check there is no expression in E. Coli sample using PA14 reference
print("Total number of genes:", len(expression_pa14_df.loc["SRR13234437"]))
print("Percent of 0 expressed genes:", 
      (expression_pa14_df.loc["SRR13234437"] == 0).sum()/len(expression_pa14_df.loc["SRR13234437"]))

print("Percent of 0 expressed genes:", 
      ((expression_pa14_df.loc["SRR13234437"] > 0) & (expression_pa14_df.loc["SRR13234437"] < 1000))
       .sum()/len(expression_pa14_df.loc["SRR13234437"]))


# In[47]:


EColi_pa14_ref_series = expression_pa14_df.loc["SRR13234437"]
nonzero_genes_EColi_pa14_ref = EColi_pa14_ref_series.iloc[EColi_pa14_ref_series.nonzero()]
sns.distplot(nonzero_genes_EColi_pa14_ref, kde=False)

# Get gene with the max expression
EColi_pa14_ref_series[EColi_pa14_ref_series > 290000]


# **Takeaway:**
# 
# * Most genes were found to have 0 expression, as expected
# * Looks like the outlier gene encodes the elongation Tu which appears to be conserved across prokaryotes: https://pubmed.ncbi.nlm.nih.gov/6985898/

# ## Positive case: PAO1
# 
# * Check that PAO1 sample (SRR13160334) aligns well to PAO1 and phage reference. 
# * Check that PAO1 sample has expression at LUZ19 (NC_010326.1 in phage genome)

# In[10]:


# Check there is expression in PAO1 sample using PAO1 reference
print("Total number of genes:", len(expression_pao1_df.loc["SRR13160334"]))
print("Percent of 0 expressed genes:", 
      (expression_pao1_df.loc["SRR13160334"] > 0).sum()/len(expression_pao1_df.loc["SRR13160334"]))


# In[58]:


# Plot the expression of PAO1 nonzero genes
pao1_sample_pao1_ref_series = expression_pao1_df.loc["SRR13160334"]
nonzero_genes_pao1_sample_pao1_ref = pao1_sample_pao1_ref_series.iloc[pao1_sample_pao1_ref_series.nonzero()]
sns.distplot(nonzero_genes_pao1_sample_pao1_ref, kde=False)

# Get gene with the max expression

len(pao1_sample_pao1_ref_series[(pao1_sample_pao1_ref_series < 1000) & pao1_sample_pao1_ref_series > 0])/len(pao1_sample_pao1_ref_series)


# In[13]:


# Check there is expression in PAO1 sample using phage reference gene: NC_010326.1
expression_phage_df.loc["SRR13160334", "NC_010326.1"]


# **Takeaway:**
# 
# * Most genes in PAO1 sample are nonzero, including phage gene that the sample is known to contain -- as expected
# * Is the TPM < 1000 expected?
# * Why is TPM of phage gene higher than other PAO1 genes?

# ## Positive case: PA14
# 
# * Check that PA14 sample (ERR3642743) aligns well to PA14 and phage references. 
# * Check that PA14 sample has expression at Pf5 (PA14_49010, PA14_49020 in PA14 genome) 
# * PA14_49010, PA14_49020 = PGD1658748, PGD1658750

# In[14]:


# Check there is expression in PA14 sample using PA14 reference
print("Total number of genes:", len(expression_pa14_df.loc["ERR3642743"]))
print("Percent of 0 expressed genes:", 
      (expression_pa14_df.loc["ERR3642743"] > 0).sum()/len(expression_pa14_df.loc["ERR3642743"]))
sns.distplot(expression_pa14_df.loc["ERR3642743"], kde=False)


# In[60]:


# Plot the expression of PAO1 nonzero genes
pa14_sample_pa14_ref_series = expression_pa14_df.loc["ERR3642743"]
nonzero_genes_pa14_sample_pa14_ref = pa14_sample_pa14_ref_series.iloc[pa14_sample_pa14_ref_series.nonzero()]
sns.distplot(nonzero_genes_pa14_sample_pa14_ref, kde=False)

# Get gene with the max expression

len(pa14_sample_pa14_ref_series[(pa14_sample_pa14_ref_series < 1000) & pa14_sample_pa14_ref_series > 0])/len(pa14_sample_pa14_ref_series)


# In[15]:


# Check there is expression in PA14 sample using PA14 reference genes: PA14_49010, PA14_49020
expression_pa14_df.loc["ERR3642743", ["PGD1658748", "PGD1658750"]]


# **Takeaways:**
# 
# * Most genes in PA14 sample are nonzero, including phage gene that the sample is known to contain -- as expected
# * Is the TPM < 1000 expected?
