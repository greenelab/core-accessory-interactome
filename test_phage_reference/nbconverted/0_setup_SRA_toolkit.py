
# coding: utf-8

# # Setup SRA toolkit
# 
# This notebook only needs to run once to download and setup SRA toolkit

# In[ ]:


import os


# In[ ]:


# Download latest version of compiled binaries of NCBI SRA toolkit 
if not os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
    get_ipython().system(' wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"')


# In[ ]:


# Extract tar.gz file 
if os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
    get_ipython().system(' tar -xzf sratoolkit.current-centos_linux64.tar.gz')

# add binaries to path using export path or editing ~/.bashrc file
get_ipython().system(' export PATH=$PATH:sratoolkit.2.10.7-centos_linux64/bin')

# Now SRA binaries added to path and ready to use

