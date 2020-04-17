#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import random
import statistics
from plotnine import (ggplot,
                      aes,
                     geom_histogram,
                     geom_density,
                     geom_col,
                     theme,
                     element_text)


# In[2]:


# Input
base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))

real_mod_file = os.path.join(
    base_dir,
    "pilot_individual_experiment",
    "data",
   "networks",
    "selected_modules.tsv")

shuffled_mod_file = os.path.join(
    base_dir,
    "pilot_individual_experiment",
    "data",
   "networks",
    "shuffled_selected_modules.tsv")

gene_annot_file = os.path.join(
    base_dir,
    "pilot_individual_experiment",
    "data",
    "annotations",
    "selected_gene_annotations.txt")


# In[3]:


# Read data
real_mods = pd.read_table(
    real_mod_file,
    header=0,
    sep='\t',
    index_col=0)

shuffled_mods = pd.read_table(
    shuffled_mod_file,
    header=0,
    sep='\t',
    index_col=0)

gene_annot = pd.read_table(
    gene_annot_file,
    header=0,
    sep='\t',
    index_col=0)

real_mods.head()


# In[4]:


shuffled_mods.head()


# In[5]:


gene_annot.head()


# ### Total number of network modules found

# In[6]:


# Total number of modules
uniq_real_mods = real_mods["dynamicColors"].unique()
uniq_shuffled_mods = shuffled_mods["dynamicColors"].unique()

if len(uniq_real_mods) == 1 and uniq_real_mods == "grey":
    tot_real_mods = 0
else:
    tot_real_mods = len(uniq_real_mods)
    
if len(uniq_shuffled_mods) == 1 and uniq_shuffled_mods == "grey":
    tot_shuffled_mods = 0
else:
    tot_shuffled_mods = len(uniq_shuffled_mods)
    
print(tot_real_mods)
print(tot_shuffled_mods)


# ### Size of modules

# In[7]:


# Number of genes per module
real_mod_groups = real_mods.groupby(["dynamicColors"]).groups

size_real_mods = []
for k,v in real_mod_groups.items():
    size_real_mods.append(len(v))
    
size_real_mods_df = pd.DataFrame(size_real_mods,
                                columns=['num_genes'])


# In[8]:


# Median size module
statistics.median(size_real_mods)


# In[9]:


# Distribution of module sizes
g = ggplot(size_real_mods_df, aes(x="num_genes"))
g += geom_density()
g


# ### Modules with both core and accessory gene interactions

# In[10]:


# For each module determine if there are any with both core and accessory genes
mods_with_core_and_acc = []
for k,v in real_mod_groups.items():
    mod_gene_annotations = gene_annot.loc[list(v)]['annotation'].unique()
    if len(mod_gene_annotations) > 1:
        mods_with_core_and_acc.append(k)
print(len(mods_with_core_and_acc))
print(mods_with_core_and_acc)


# In[11]:


# Merge annotation and module information per gene
gene_annot_mods = gene_annot.merge(real_mods, left_index=True, right_index=True)

gene_annot_mods.head()


# In[12]:


# Create df 
mod_names = []
core_acc_names = []
val = []
for k,v in real_mod_groups.items():
    for gene_cat in ["core", "accessory"]:
        mod_names.append(k)
        core_acc_names.append(gene_cat)
        lst_gene_annotations = list(gene_annot.loc[list(v)]['annotation'])
        val.append(lst_gene_annotations.count(gene_cat))


# In[13]:


df = pd.DataFrame(data = list(zip(mod_names, core_acc_names, val)),
                  columns=["modules", "gene annot", "value"])

df.head()


# In[14]:


# Modules with only core genes
df[df['value'] == 0]


# In[15]:


(ggplot(df, aes(x='modules', y='value', fill='gene annot'))
 + geom_col()
 + theme(axis_text_x=element_text(rotation=90, hjust=1))
)


# In[16]:


# Next questions to ask
# Pick modules with only core
# Pick modules with mixed core and accessory
# What types of genes are in these two modules?
# Enrichment of operons, pathways, can we get a set of operons/pathways known to have both core and accessory? 


# ### Manually explore genes in modules
# 
# Cherrry pick core-only module and core-accessory module and search for literature about what is known about genes

# In[21]:


# Gene number to gene name file
gene_name_file = os.path.join(
    base_dir,
    "pilot_individual_experiment",
    "data",
    "annotations",
    "Pseudomonas_aeruginosa_PAO1_107.csv")


# In[22]:


# Read gene number to name mapping
gene_name_mapping = pd.read_table(
    gene_name_file,
    header=0,
    sep=',',
    index_col=0)

gene_name_mapping = gene_name_mapping[["Locus Tag", "Name"]]

gene_name_mapping.set_index("Locus Tag", inplace=True)
gene_name_mapping.head()


# In[23]:


# Format gene numbers to remove extraneous quotes
gene_number = gene_name_mapping.index
gene_name_mapping.index = gene_number.str.strip("\"")

gene_name_mapping.head()


# In[24]:


# Map gene numbers to names
def get_gene_names(gene_id_list):    
    gene_names = []
    for gene_id in gene_id_list:
        gene_name = gene_name_mapping.loc[gene_id]
        if gene_name.isnull()[0]:
            # If gene name does not exist
            # Use gene number
            gene_names.append(gene_id)
        else:
            gene_names.append(gene_name[0])
    return gene_names


# In[43]:


# Select core-only module
#core_only_gene_ids = list(real_mod_groups['darkgreen'])
core_only_gene_ids = list(real_mod_groups['yellowgreen'])

core_only_annot = gene_annot.loc[core_only_gene_ids]['annotation']
core_only_gene_names = get_gene_names(core_only_gene_ids)

core_only_df = pd.DataFrame(data = list(zip(core_only_gene_names, core_only_annot)),
                           columns=["gene names", "annot"],
                          index = core_only_gene_ids)

core_only_df.head(100)


# In[40]:


# Select core-oacc module
core_acc_gene_ids = list(real_mod_groups['grey60'])
#core_acc_gene_ids = list(real_mod_groups['lightyellow'])

core_acc_annot = gene_annot.loc[core_acc_gene_ids]['annotation']
core_acc_gene_names = get_gene_names(core_acc_gene_ids)

core_acc_df = pd.DataFrame(data = list(zip(core_acc_gene_names, core_acc_annot)),
                           columns=["gene names", "annot"],
                          index = core_acc_gene_ids)

core_acc_df.head(30)

