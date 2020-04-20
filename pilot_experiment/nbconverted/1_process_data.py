#!/usr/bin/env python
# coding: utf-8

# # Generate gene-gene network modules
# This notebook generate gene-gene network modules for an exploratory analysis.
# 
# The objective is____

# In[1]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
import pandas as pd
import os
import argparse
from functions import process_data


# In[2]:


base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))


# In[3]:


# Input files
normalized_data_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
   "input",
    "train_set_normalized.pcl")

metadata_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "annotations",
    "sample_annotations.tsv")

# Load in annotation file
# Annotation file contains the list of all PAO1 specific genes
gene_mapping_file = os.path.join(
    base_dir,
    "pilot_experiment",
    "data",
    "annotations",
    "PAO1_ID_PA14_ID.csv")


# In[4]:


# Select specific experiment
lst_experiments = ["E-GEOD-8083",
                   "E-GEOD-29789",
                   "E-GEOD-48982",
                   "E-GEOD-24038",
                   "E-GEOD-29879",
                   "E-GEOD-49759"]


# In[5]:


# Output files
selected_data_file = os.path.join(
        base_dir,
        "pilot_experiment",
        "data",
        "input",
        "selected_normalized_data.tsv")

shuffled_selected_data_file = os.path.join(
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


# # Select subset of samples
# 
# Select experiments that contain either PAO1 or PA14 strains. We will use only PAO1 and PA14 strains as a first pass because these two strains are the most common and well studied *P. aeruginosa* strains, therefore we will be able to verify the resulting gene-gene interactions with those found in the literature.

# In[6]:


process_data.select_expression_data(normalized_data_file,
                                   metadata_file,
                                   lst_experiments,
                                   selected_data_file)


# # Permute selected expression data
# This permuted version will serve as a baseline for our analysis

# In[7]:


process_data.permute_expression_data(selected_data_file,
                                     shuffled_selected_data_file)


# # Annotate genes as core and accessory
# 
# Annotate genes as either **core** if PAO1 gene is homologous to PA14 gene, or **accessory** if there does not exist a homolog. 
# 
# These homologous mappings are based on the [Bactome database](https://bactome.helmholtz-hzi.de/cgi-bin/h-pange.cgi?STAT=1&Gene=PA0135)

# In[8]:


process_data.annotate_genes(selected_data_file,
                            gene_mapping_file,
                            gene_annot_file)


# # Network construction and module detection
# 
# Networks provide a straightforward representation of interactions between the nodes. A node corresponds to the gene expression profile of a given gene. Nodes are connected if they have a significant pairwise expression profile association across the environmental perturbations (cell- or tissue- samples). It is standard to use the (Pearson) correlation coefficient as a co-expression measure, e.g., the absolute value of Pearson correlation. 

# ## Get parameters for network generation
# 
# **Question:** How do we pick a threshold to determine what Pearson correlation score is sufficient to say that 2 nodes are associated?
# 
# A 'hard threshold' may lead to loss of information and sensitity.
# Instead, a 'soft thresholding' is proposed. Soft thresholding weighs each connection of a float [0,1]   
# 
# **Reference:**
# * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.471.9599&rep=rep1&type=pdf

# In[9]:


get_ipython().run_cell_magic('R', '-i selected_data_file', '# Get threshold param for true gene expression data\nsource("functions/network_utils.R")\n\nget_threshold(selected_data_file)')


# In[10]:


get_ipython().run_cell_magic('R', '-i shuffled_selected_data_file', '# Get threshold param for true gene expression data\nsource("functions/network_utils.R")\n\nget_threshold(shuffled_selected_data_file)')


# In[11]:


# Prompt user to choose threshold params
power_param_true = int(input("Treshold for true data:"))


# In[12]:


# Prompt user to choose threshold params
power_param_shuffled = int(input("Threshold for permuted data:"))


# ## Generate network modules

# In[13]:


# Output module files
gene_modules_file = os.path.join(
        base_dir,
        "pilot_experiment",
        "data",
        "networks",
        "selected_modules.tsv")

shuffled_gene_modules_file = os.path.join(
        base_dir,
        "pilot_experiment",
        "data",
        "networks",
        "shuffled_selected_modules.tsv")


# In[14]:


get_ipython().run_cell_magic('R', '-i power_param_true -i selected_data_file -i gene_modules_file', '# Generate network modules using threshold params selected for true expression data\nsource("functions/network_utils.R")\n\ngenerate_network_modules(power_param_true,\n                         selected_data_file,\n                         gene_modules_file)')


# In[15]:


get_ipython().run_cell_magic('R', '-i power_param_shuffled -i shuffled_selected_data_file -i shuffled_gene_modules_file', '# Generate network modules using threshold params selected for permuted expression data\nsource("functions/network_utils.R")\n\ngenerate_network_modules(power_param_shuffled,\n                         shuffled_selected_data_file,\n                         shuffled_gene_modules_file)')

