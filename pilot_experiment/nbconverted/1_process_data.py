
# coding: utf-8

# # Process data
# This notebook processes the expression data that will be used in this pilot analysis. Specifically this notebook performs the following steps:
# 
# 1. Selects a small subset of expression data and outputs this dataset to file
# 2. Permutes the subsetted data to use as a control and outputs this dataset to file
# 3. Generates a mapping between *P. aeruginosa* gene id (PA####) and core, accessory label

# In[1]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
import pandas as pd
import os
import argparse
import numpy as np
from functions import process_data

np.random.seed(123)


# In[2]:


base_dir = os.path.abspath(os.path.join(os.getcwd(),"../"))


# ### About input data
# Normalized expression data is the *P. aeruginosa* compendium from [Tan et. al.](https://msystems.asm.org/content/1/1/e00025-15). The dataset can be found in the associated [ADAGE github repository](https://github.com/greenelab/adage/blob/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl).
# 
# The corresponding metadata was downloaded from the [ADAGE website](https://adage.greenelab.com/#/download).

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


# Select specific experiments
# In this case we selected 6 experiments (3 experiments PAO1 strains and the other 3 experiments contain PA14 strains).
# We will use only PAO1 and PA14 strains as a first pass because these two strains are the most common and well studied
# P. aeruginosa strains and therefore we will be able to verify the resulting gene-gene interactions with those found
# in the literature.
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
# Select subset of experiments to use

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

