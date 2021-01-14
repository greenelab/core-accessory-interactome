#!/usr/bin/env python
# coding: utf-8

# # Setup SRA toolkit
# 
# This notebook only needs to run once to download and setup SRA toolkit, which we use to download transcriptome data from SRA to process with Salmon

# In[1]:


import os


# In[2]:


# Download latest version of compiled binaries of NCBI SRA toolkit 
get_ipython().system(' wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz')


# In[4]:


# Extract tar.gz file 
get_ipython().system('tar -vxzf sratoolkit.tar.gz')


# In[22]:


# add binaries to path using export path or editing ~/.bashrc file
# Needs to be run in terminal for some reason
get_ipython().system('export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin')

# Now SRA binaries added to path and ready to use

