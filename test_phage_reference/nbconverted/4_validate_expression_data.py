#!/usr/bin/env python
# coding: utf-8

# # Validate expression data
# 
# This notebook is checking that the resulting expression data is as expected using a few test samples.

# In[1]:


import seaborn as sns
import pandas as pd
import numpy as np
from core_acc_modules import paths


# In[2]:


# Load expression data
expression_pao1_phage_df = pd.read_csv(paths.PAO1_PHAGE_GE, sep='\t', index_col=0, header=0)
expression_pa14_phage_df = pd.read_csv(paths.PA14_PHAGE_GE, sep='\t', index_col=0, header=0)


# In[3]:


expression_pao1_phage_df.head()


# In[4]:


expression_pa14_phage_df.head()


# ## Some checks

# In[5]:


# Are there any genes that are 0 across all samples?
expression_pao1_phage_df.loc[:,expression_pao1_phage_df.sum()==0] 


# In[6]:


expression_pa14_phage_df.loc[:,expression_pa14_phage_df.sum()==0] 


# In[7]:


# What does the distribution of activity look like? Are there outliers?
#sns.distplot(np.log10(expression_pao1_phage_df.mean()))


# In[8]:


#sns.distplot(expression_pa14_phage_df.mean(), kde=False)


# ## Negative case
# 
# Check that *E. Coli* sample (SRR13234437) poorly aligns with all *P. aeruginosa* reference genomes including PAO1+phage, PA14+phage.

# In[9]:


# Check there is no expression in E. Coli sample using phage reference
print("Total number of genes:", len(expression_pao1_phage_df.loc["SRR13234437"]))
print("Percent of 0 expressed genes:",
      (expression_pao1_phage_df.loc["SRR13234437"] == 0).sum()/len(expression_pao1_phage_df.loc["SRR13234437"]))


# In[10]:


EColi_pao1_phage_ref_series = expression_pao1_phage_df.loc["SRR13234437"]
nonzero_genes_EColi_pao1_phage_ref = EColi_pao1_phage_ref_series.iloc[EColi_pao1_phage_ref_series.nonzero()]

# Plot
fig2 = sns.distplot(nonzero_genes_EColi_pao1_phage_ref, kde=False)
fig2.set_title("Expression of E. Coli sample mapped to PAO1+phage reference (TPM)")

# Get gene with the max expression
# What is the outlier one?


# In[11]:


# If we remove the outlier gene, what does the expression look like
# Does this distribution look correct?
# What are the genes outside of this distribution? Why are they outside
threshold = 1000
majority_genes_EColi_pao1_phage_ref = nonzero_genes_EColi_pao1_phage_ref[nonzero_genes_EColi_pao1_phage_ref < threshold]

print(f"Percent of nonzero genes below {threshold}:", 
      (len(majority_genes_EColi_pao1_phage_ref))/len(nonzero_genes_EColi_pao1_phage_ref))

fig5 = sns.distplot(majority_genes_EColi_pao1_phage_ref, kde=False)


# In[12]:


# Check there is no expression in E. Coli sample using PAO1 reference
print("Total number of genes:", len(expression_pa14_phage_df.loc["SRR13234437"]))
print("Percent of 0 expressed genes:", 
      (expression_pa14_phage_df.loc["SRR13234437"] == 0).sum()/len(expression_pa14_phage_df.loc["SRR13234437"]))


# In[13]:


EColi_pa14_phage_ref_series = expression_pa14_phage_df.loc["SRR13234437"]
nonzero_genes_EColi_pa14_phage_ref = EColi_pa14_phage_ref_series.iloc[EColi_pa14_phage_ref_series.nonzero()]

# Plot
fig2 = sns.distplot(nonzero_genes_EColi_pa14_phage_ref, kde=False)
fig2.set_title("Expression of E. Coli sample mapped to PA14+phage reference (TPM)")

# Get gene with the max expression
# What is the outlier one?


# In[14]:


# If we remove the outlier gene, what does the expression look like
# Does this distribution look correct?
# What are the genes outside of this distribution? Why are they outside
threshold = 1000
majority_genes_EColi_pa14_phage_ref = nonzero_genes_EColi_pa14_phage_ref[nonzero_genes_EColi_pa14_phage_ref < threshold]

print(f"Percent of nonzero genes below {threshold}:", 
      (len(majority_genes_EColi_pa14_phage_ref))/len(nonzero_genes_EColi_pa14_phage_ref))

fig5 = sns.distplot(majority_genes_EColi_pa14_phage_ref, kde=False)


# In[15]:


# Get outlier gene
print(nonzero_genes_EColi_pao1_phage_ref[nonzero_genes_EColi_pao1_phage_ref == nonzero_genes_EColi_pao1_phage_ref.max()])
print(nonzero_genes_EColi_pa14_phage_ref[nonzero_genes_EColi_pa14_phage_ref == nonzero_genes_EColi_pa14_phage_ref.max()])


# ## Positive case: PAO1
# 
# Given a PAO1 sample (SRR13160334), we want to check that: 
# * It aligns well to PAO1 reference. 
# * It has expression at LUZ19 phage (NC_010326.1 in phage genome), since this sample was known to contain this phage gene.

# In[16]:


# Check there is expression in PAO1 sample using PAO1 reference
print("Total number of genes:", len(expression_pao1_phage_df.loc["SRR13160334"]))
print("Percent of positive-expressed genes:", 
      (expression_pao1_phage_df.loc["SRR13160334"] > 0).sum()/len(expression_pao1_phage_df.loc["SRR13160334"]))


# In[17]:


# Plot the expression of PAO1 nonzero genes
pao1_sample_pao1_phage_ref_series = expression_pao1_phage_df.loc["SRR13160334"]
nonzero_genes_pao1_sample_pao1_phage_ref = pao1_sample_pao1_phage_ref_series.iloc[pao1_sample_pao1_phage_ref_series.nonzero()]

# Plot
fig3 = sns.distplot(nonzero_genes_pao1_sample_pao1_phage_ref, kde=False)
fig3.set_title("Expression of PAO1 sample mapped to PAO1+phage reference (TPM)")


# In[18]:


# If we remove the outlier gene, what does the expression look like
# Does this distribution look correct?
# What are the genes outside of this distribution? Why are they outside
threshold = 1000
majority_genes_pao1_sample_pao1_phage_ref = nonzero_genes_pao1_sample_pao1_phage_ref[nonzero_genes_pao1_sample_pao1_phage_ref < threshold]

print(f"Percent of nonzero genes below {threshold}:", 
      (len(majority_genes_pao1_sample_pao1_phage_ref))/len(nonzero_genes_pao1_sample_pao1_phage_ref))

fig5 = sns.distplot(nonzero_genes_pao1_sample_pao1_phage_ref[nonzero_genes_pao1_sample_pao1_phage_ref < threshold], kde=False)


# In[19]:


# Check there is expression in PAO1 sample using phage reference gene: NC_010326.1
expression_pao1_phage_df.loc["SRR13160334", "NC_010326.1"]


# In[20]:


# Look at outlier genes
nonzero_genes_pao1_sample_pao1_phage_ref[nonzero_genes_pao1_sample_pao1_phage_ref > threshold]


# ## Positive case: PA14
# 
# Given a PA14 sample (ERR3642743), we want to check that:
# * It aligns well to PA14 reference. 
# * It has expression at Pf5 (PA14_49010, PA14_49020 in PA14 genome) which the sample is known to contain
# 
# Implementation note:
# * Here is the mapping from PA14 id to sequence id: PA14_49010, PA14_49020 = PGD1658748, PGD1658750

# In[21]:


# Check there is expression in PA14 sample using PA14 reference
print("Total number of genes:", len(expression_pa14_phage_df.loc["ERR3642743"]))
print("Percent of positive-expressed genes:", 
      (expression_pa14_phage_df.loc["ERR3642743"] > 0).sum()/len(expression_pa14_phage_df.loc["ERR3642743"]))


# In[22]:


# Plot the expression of PAO1 nonzero genes
pa14_sample_pa14_phage_ref_series = expression_pa14_phage_df.loc["ERR3642743"]
nonzero_genes_pa14_sample_pa14_phage_ref = pa14_sample_pa14_phage_ref_series.iloc[pa14_sample_pa14_phage_ref_series.nonzero()]

# Plot
fig4 = sns.distplot(nonzero_genes_pa14_sample_pa14_phage_ref, kde=False)
fig4.set_title("Expression of PA14 sample mapped to PA14+phage reference (TPM)")

# Get gene with the max expression
#len(pa14_sample_pa14_ref_series[(pa14_sample_pa14_ref_series < 1000) & pa14_sample_pa14_ref_series > 0])/len(pa14_sample_pa14_ref_series)


# In[23]:


# If we remove the outlier gene, what does the expression look like
# Does this distribution look correct?
# What are the genes outside of this distribution? Why are they outside
threshold = 1000
majority_genes_pa14_sample_pa14_phage_ref = nonzero_genes_pa14_sample_pa14_phage_ref[nonzero_genes_pa14_sample_pa14_phage_ref < threshold]

print(f"Percent of nonzero genes below {threshold}:", 
      (len(majority_genes_pa14_sample_pa14_phage_ref))/len(nonzero_genes_pa14_sample_pa14_phage_ref))

fig5 = sns.distplot(nonzero_genes_pa14_sample_pa14_phage_ref[nonzero_genes_pa14_sample_pa14_phage_ref < threshold], kde=False)


# In[24]:


# Check there is expression in PA14 sample using PA14 reference genes: PA14_49010, PA14_49020
expression_pa14_phage_df.loc["ERR3642743", ["PGD1658748", "PGD1658750"]]


# In[25]:


# Look at outlier genes
nonzero_genes_pa14_sample_pa14_phage_ref[nonzero_genes_pa14_sample_pa14_phage_ref > threshold]


# **Takeaways:**
# 
# As expected,
# * 92%, 96% of genes in E. Coli samples had 0 expression using PAO1+phage reference and PA14+phage reference respectively
# * 78% of genes in PAO1 sample had positive expression using PAO1+phage reference, including phage gene it was known to contain
# * 83% of genes in PA14 sample had positive expression using PA14+phage reference, including phage gene it was known to contain
# * Of the genes that had nonzero expression, most (97-98%) were < 1000 TPM nonzero genes were found. Genes exceeding this threshold are rRNAs, tRNAs
# 
# On the other hand we found some unexpected mapping_rates.
# * Some rates are nearly equal when mapping a PAO1 sample using the PAO1 reference versus PA14 reference
# * Overall many mappings are close to around 50% which seems low
# * No clear trend relating processing methods to the samples that had unexpected mapping patterns
# 
# (Slide 10): https://docs.google.com/presentation/d/1QqyNcvyLF83bdxgqyGFV2fLSfSs2yP-9eg8fQeSN5mU/edit?usp=sharing
# 
# TO DO:
# * BLAST sequences that should differ between PAO1 and PA14 against PAO1 and PA14 references to make sure that the PA14 reference actually contains PA14 sequences.
