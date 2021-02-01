
# coding: utf-8

# # Clustering by metadata
# 
# This notebook examines how the data clusters by metadata. We would expect that the data cluster by strain type (PAO1 and PA14).
# 
# Other metadata?

# In[1]:


import os
import pandas as pd
from sklearn.decomposition import PCA
import umap
import plotnine as pn
from core_acc_modules import paths_corr


# In[2]:


# Load expression data
expression_pao1_df = pd.read_csv(paths_corr.PAO1_GE, sep='\t', index_col=0, header=0)
expression_pa14_df = pd.read_csv(paths_corr.PA14_GE, sep='\t', index_col=0, header=0)


# In[3]:


# Load metadata
pao1_metadata = pd.read_csv(paths_corr.PAO1_METADATA, index_col=0, header=0)
pa14_metadata = pd.read_csv(paths_corr.PA14_METADATA, index_col=0, header=0)


# In[4]:


# Select meatadata variable to color by
#metadata_variable = 'genotype'
metadata_variable = 'processing'


# In[5]:


print(expression_pao1_df.shape)
expression_pao1_df.head()


# In[6]:


print(expression_pa14_df.shape)
expression_pa14_df.head()


# ## Parse metadata

# ## Plot in PCA space

# In[7]:


# Embed expression data into low dimensional space
pca = PCA(n_components=2)
model = pca.fit(expression_pao1_df)
pao1_encoded = model.transform(expression_pao1_df)

pao1_encoded_df = pd.DataFrame(data=pao1_encoded,
                               index=expression_pao1_df.index,
                               columns=['1','2'])

# Add label
if metadata_variable == "genotype":
    pa14_ids = list(pa14_metadata.index)
    pao1_encoded_df['genotype'] = 'PAO1'
    pao1_encoded_df.loc[pa14_ids,'genotype'] = 'PA14'
elif metadata_variable == "processing":
    pao1_ids = list(pao1_metadata.index)
    pa14_ids = list(pa14_metadata.index)
    pao1_encoded_df.loc[pao1_ids,"processing"] = pao1_metadata.loc[pao1_ids,"LibrarySelection"].values
    pao1_encoded_df.loc[pa14_ids,"processing"] = pa14_metadata.loc[pa14_ids,"LibrarySelection"].values

pao1_encoded_df.head()


# In[8]:


model = pca.fit(expression_pa14_df)
pa14_encoded = model.transform(expression_pa14_df)

pa14_encoded_df = pd.DataFrame(data=pa14_encoded,
                               index=expression_pa14_df.index,
                               columns=['1','2'])

# Add label
if metadata_variable == "genotype":
    pa14_ids = list(pa14_metadata.index)
    pa14_encoded_df['genotype'] = 'PAO1'
    pa14_encoded_df.loc[pa14_ids,'genotype'] = 'PA14'
elif metadata_variable == "processing":
    pao1_ids = list(pao1_metadata.index)
    pa14_ids = list(pa14_metadata.index)
    pa14_encoded_df.loc[pao1_ids,"processing"] = pao1_metadata.loc[pao1_ids,"LibrarySelection"].values
    pa14_encoded_df.loc[pa14_ids,"processing"] = pa14_metadata.loc[pa14_ids,"LibrarySelection"].values

pa14_encoded_df.head()


# In[9]:


# Plot PAO1
fig = pn.ggplot(pao1_encoded_df, pn.aes(x='1', y='2'))
fig += pn.geom_point(pn.aes(color=metadata_variable), alpha=0.5)
fig += pn.labs(x ='PC 1',
            y = 'PC 2',
            title = 'RNA-seq expression using PAO1 reference')
fig += pn.theme_bw()
fig += pn.theme(
    legend_title_align = "center",
    plot_background=pn.element_rect(fill='white'),
    legend_key=pn.element_rect(fill='white', colour='white'), 
    legend_title=pn.element_text(family='sans-serif', size=15),
    legend_text=pn.element_text(family='sans-serif', size=12),
    plot_title=pn.element_text(family='sans-serif', size=15),
    axis_text=pn.element_text(family='sans-serif', size=12),
    axis_title=pn.element_text(family='sans-serif', size=15)
    )
fig += pn.guides(colour=pn.guide_legend(override_aes={'alpha': 1}))

print(fig)


# In[10]:


# Plot PA14
fig = pn.ggplot(pa14_encoded_df, pn.aes(x='1', y='2'))
fig += pn.geom_point(pn.aes(color=metadata_variable), alpha=0.5)
fig += pn.labs(x ='PC 1',
            y = 'PC 2',
            title = 'RNA-seq expression using PA14 reference')
fig += pn.theme_bw()
fig += pn.theme(
    legend_title_align = "center",
    plot_background=pn.element_rect(fill='white'),
    legend_key=pn.element_rect(fill='white', colour='white'), 
    legend_title=pn.element_text(family='sans-serif', size=15),
    legend_text=pn.element_text(family='sans-serif', size=12),
    plot_title=pn.element_text(family='sans-serif', size=15),
    axis_text=pn.element_text(family='sans-serif', size=12),
    axis_title=pn.element_text(family='sans-serif', size=15)
    )
fig += pn.guides(colour=pn.guide_legend(override_aes={'alpha': 1}))

print(fig)


# ## Plot in UMAP

# In[11]:


# Embed expression data into low dimensional space
model = umap.UMAP(random_state=123).fit(expression_pao1_df)
pao1_encoded = model.transform(expression_pao1_df)

pao1_encoded_df = pd.DataFrame(data=pao1_encoded,
                               index=expression_pao1_df.index,
                               columns=['1','2'])

# Add label
if metadata_variable == "genotype":
    pa14_ids = list(pa14_metadata.index)
    pao1_encoded_df['genotype'] = 'PAO1'
    pao1_encoded_df.loc[pa14_ids,'genotype'] = 'PA14'
elif metadata_variable == "processing":
    pao1_ids = list(pao1_metadata.index)
    pa14_ids = list(pa14_metadata.index)
    pao1_encoded_df.loc[pao1_ids,"processing"] = pao1_metadata.loc[pao1_ids,"LibrarySelection"].values
    pao1_encoded_df.loc[pa14_ids,"processing"] = pa14_metadata.loc[pa14_ids,"LibrarySelection"].values

pao1_encoded_df.head()


# In[12]:


# Embed expression data into low dimensional space
model = umap.UMAP(random_state=123).fit(expression_pa14_df)
pa14_encoded = model.transform(expression_pa14_df)

pa14_encoded_df = pd.DataFrame(data=pa14_encoded,
                               index=expression_pa14_df.index,
                               columns=['1','2'])

# Add label
if metadata_variable == "genotype":
    pa14_ids = list(pa14_metadata.index)
    pa14_encoded_df['genotype'] = 'PAO1'
    pa14_encoded_df.loc[pa14_ids,'genotype'] = 'PA14'
elif metadata_variable == "processing":
    pao1_ids = list(pao1_metadata.index)
    pa14_ids = list(pa14_metadata.index)
    pa14_encoded_df.loc[pao1_ids,"processing"] = pao1_metadata.loc[pao1_ids,"LibrarySelection"].values
    pa14_encoded_df.loc[pa14_ids,"processing"] = pa14_metadata.loc[pa14_ids,"LibrarySelection"].values

pa14_encoded_df.head()


# In[13]:


# Plot PAO1
fig = pn.ggplot(pao1_encoded_df, pn.aes(x='1', y='2'))
fig += pn.geom_point(pn.aes(color=metadata_variable), alpha=0.5)
fig += pn.labs(x ='UMAP 1',
            y = 'UMAP 2',
            title = 'RNA-seq expression using PAO1 reference')
fig += pn.theme_bw()
fig += pn.theme(
    legend_title_align = "center",
    plot_background=pn.element_rect(fill='white'),
    legend_key=pn.element_rect(fill='white', colour='white'), 
    legend_title=pn.element_text(family='sans-serif', size=15),
    legend_text=pn.element_text(family='sans-serif', size=12),
    plot_title=pn.element_text(family='sans-serif', size=15),
    axis_text=pn.element_text(family='sans-serif', size=12),
    axis_title=pn.element_text(family='sans-serif', size=15)
    )
fig += pn.guides(colour=pn.guide_legend(override_aes={'alpha': 1}))

print(fig)


# In[14]:


# Plot PA14
fig = pn.ggplot(pa14_encoded_df, pn.aes(x='1', y='2'))
fig += pn.geom_point(pn.aes(color=metadata_variable), alpha=0.5)
fig += pn.labs(x ='UMAP 1',
            y = 'UMAP 2',
            title = 'RNA-seq expression using PA14 reference')
fig += pn.theme_bw()
fig += pn.theme(
    legend_title_align = "center",
    plot_background=pn.element_rect(fill='white'),
    legend_key=pn.element_rect(fill='white', colour='white'), 
    legend_title=pn.element_text(family='sans-serif', size=15),
    legend_text=pn.element_text(family='sans-serif', size=12),
    plot_title=pn.element_text(family='sans-serif', size=15),
    axis_text=pn.element_text(family='sans-serif', size=12),
    axis_title=pn.element_text(family='sans-serif', size=15)
    )
fig += pn.guides(colour=pn.guide_legend(override_aes={'alpha': 1}))

print(fig)


# **Takeaway:**
# 
# * Looks like there is a pretty clear separation between PAO1 and PA14 strains using PAO1 reference. Looks like library selection could explain some of the variance found in using the PA14 reference. 
# 
# * We will need to determine if this is something we will need to correct in our analysis or if these sources of variance will become noise with a larger dataset (i.e. the compendium that Georgia processed)
