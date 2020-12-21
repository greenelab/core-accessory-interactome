
# coding: utf-8

# # Align and quantify samples
# 
# This notebook aligns test samples against the phage and PAO1 reference genomes. Our goal is to test our phage reference genome alignment before we port it to the Discovery (Dartmouth computing cluster). We want to check that we are getting more expression of phage genes in a phage sample compared to a non-pseudomonas samples.

# In[17]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import shutil
import pandas as pd
import numpy as np
from core_acc_modules import paths

np.random.seed(123)


# ### Download SRA data
# 
# Note: Need to delete `sra` folder between runs otherwise `fastq-dump` will be called on all files in `sra` folder which can include more than your sra accessions.

# In[4]:


shutil.rmtree(paths.SRA_DIR)


# In[5]:


# Download sra data files
get_ipython().system(' prefetch --option-file $paths.SRA_ACC')


# ### Get FASTQ files associated with SRA downloads
# 
# The fastq files store the RNA-seq results, including: sequencing and quality scores for each base call.
# 
# Here is a nice blog to explain how to read fastq files: https://thesequencingcenter.com/knowledge-base/fastq-files/
# 
# The fastq files gives the sequence of a read at a given location. Our goal is to map these reads to a reference genome so that we can quantify the number of reads that are at a given location, to determine the level of expression.

# In[6]:


if not os.path.exists(paths.FASTQ_DIR):
    os.makedirs(paths.FASTQ_DIR)


# In[7]:


get_ipython().system('fastq-dump $paths.SRA_DIR/* --split-files --outdir $paths.FASTQ_DIR/')


# In[8]:


# Copied from https://github.com/hoganlab-dartmouth/sraProcessingPipeline/blob/5974e040c85724a8d385e53153b7707ae7c9c255/DiscoveryScripts/quantifier.py#L83

#!fastq-dump $paths_phage.SRA_DIR/* --skip-technical --readids --split-3 --clip --outdir $paths_phage.FASTQ_DIR/


# ### Quantify gene expression
# Now that we have our index built and all of our data downloaded, we’re ready to quantify our samples
# 
# **Input:**
# * Index of reference transcriptome
# * FASTQ of experimental samples
# 
# **Output:**
# 
# After the salmon commands finish running, you should have a directory named quants, which will have a sub-directory for each sample. These sub-directories contain the quantification results of salmon, as well as a lot of other information salmon records about the sample and the run. 
# 
# The main output file (called `quant.sf`). Inside the quantification file for sample DRR016125 in quants/DRR016125/quant.sf, you’ll see a TSV format file listing the name (`Name`) of each transcript, its length (`Length`), effective length (`EffectiveLength`), and its abundance in terms of Transcripts Per Million (`TPM`) and estimated number of reads (`NumReads`) originating from this transcript.
# 
# **For each sample we have read counts per gene (where the genes are based on the reference gene file provided above).** 

# #### Get quants using PAO1 reference

# In[ ]:


if not os.path.exists(paths.PAO1_QUANT):
    os.makedirs(paths.PAO1_QUANT)


# In[13]:


get_ipython().run_cell_magic('bash', '-s $paths.PAO1_QUANT $paths.FASTQ_DIR $paths.PAO1_INDEX', '\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${sample_name}"\n\nsalmon quant -i $3 -l A \\\n            -1 ${base_name}/${sample_name}_1.fastq \\\n            -2 ${base_name}/${sample_name}_2.fastq \\\n            -p 8 --validateMappings -o $1/${sample_name}_quant\ndone')


# #### Get quants using phage reference

# In[ ]:


if not os.path.exists(paths.PHAGE_QUANT):
    os.makedirs(paths.PHAGE_QUANT)


# In[12]:


get_ipython().run_cell_magic('bash', '-s $paths.PHAGE_QUANT $paths.FASTQ_DIR $paths.PHAGE_INDEX', '\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${sample_name}"\n\nsalmon quant -i $3 -l A \\\n            -1 ${base_name}/${sample_name}_1.fastq \\\n            -2 ${base_name}/${sample_name}_2.fastq \\\n            -p 8 --validateMappings -o $1/${sample_name}_quant\ndone')


# ### Consolidate sample quantification to gene expression dataframe

# In[14]:


# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PAO1_QUANT

expression_pao1_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_pao1_df.head()


# In[15]:


# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PHAGE_QUANT

expression_phage_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_phage_df.head()


# In[19]:


# Save gene expression data
expression_pao1_df.to_csv(paths.PAO1_GE, sep='\t')
expression_phage_df.to_csv(paths.PHAGE_GE, sep='\t')

