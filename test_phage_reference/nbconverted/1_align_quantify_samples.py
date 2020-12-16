#!/usr/bin/env python
# coding: utf-8

# # Aign and quantify samples
# 
# This notebook aligns test samples against the phage reference genome
# 
# *Positive test cases:*
# * ???
# 
# *Negative test cases:*
# * E. Coli sample
# * Pseudomonas sample containing only core genes

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import pandas as pd
import numpy as np
from core_acc_modules import paths_phage

np.random.seed(123)


# ### Setup SRA toolkit -- only needs to be run once

# In[2]:


# Download latest version of compiled binaries of NCBI SRA toolkit 
#if not os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
#    ! wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"


# In[3]:


# Extract tar.gz file 
#if os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
#    ! tar -xzf sratoolkit.current-centos_linux64.tar.gz

# add binaries to path using export path or editing ~/.bashrc file
#! export PATH=$PATH:sratoolkit.2.10.7-centos_linux64/bin

# Now SRA binaries added to path and ready to use


# ### Download SRA data

# In[4]:


"""# Download sra data files
! prefetch --option-file $paths.SRA_ACC """


# ### Get FASTQ files associated with SRA downloads
# 
# The fastq files store the RNA-seq results, including: sequencing and quality scores for each base call.
# 
# Here is a nice blog to explain how to read fastq files: https://thesequencingcenter.com/knowledge-base/fastq-files/
# 
# The fastq files gives the sequence of a read at a given location. Our goal is to map these reads to a reference genome so that we can quantify the number of reads that are at a given location, to determine the level of expression.

# In[5]:


#!mkdir $paths.FASTQ_DIR


# In[6]:


#!fastq-dump $paths.SRA_DIR/* --split-files --outdir $paths.FASTQ_DIR/


# ### Obtain a transcriptome and build an index
# 
# Here we are using [Salmon](https://combine-lab.github.io/salmon/)
# 
# **Input:**
# * Target transcriptome
# * This transcriptome is given to Salmon in the form of a (possibly compressed) multi-FASTA file, with each entry providing the sequence of a transcript
# * We downloaded the phage GENOMES from NCBI GenBank
# 
# **Note:** For prokaryotes, transcripts and genes have a more 1-1 mapping so we're using genes for our reference transcriptome and so we don't need to use tximport to map transcript quants to genes. 
# 
# **Output:**
# * The index is a structure that salmon uses to quasi-map RNA-seq reads during quantification
# * [Quasi-map](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985) is a way to map sequenced fragments (single or paired-end reads) to a target transcriptome. Quasi-mapping produces what we refer to as fragment mapping information. In particular, it provides, for each query (fragment), the reference sequences (transcripts), strand and position from which the query may have likely originated. In many cases, this mapping information is sufficient for downstream analysis like quantification.
# 
# *Algorithm:*
# 
# For a query read r through repeated application of: 
# 1. Determining the next hash table k-mer that starts past the current query position
# 2. Computing the maximum mappable prefix (MMP) of the query beginning with this k-mer
# 3. Determining the next informative position (NIP) by performing a longest common prefix (LCP) query on two specifically chosen suffixes in the SA

# In[7]:


#! salmon index -t $paths_phage.PAO1_REF -i $paths_phage.PAO1_INDEX


# In[8]:


get_ipython().system(' salmon index -t $paths_phage.PHAGE_REF -i $paths_phage.PHAGE_INDEX')


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

# #### Get quants using phage reference

# In[9]:


get_ipython().run_cell_magic('bash', '-s $paths_phage.PHAGE_QUANT $paths_phage.FASTQ_DIR $paths_phage.PHAGE_INDEX', 'mkdir $1\n\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${sample_name}"\n\nsalmon quant -i $3 -l A \\\n            -1 ${base_name}/${sample_name}_1.fastq \\\n            -2 ${base_name}/${sample_name}_2.fastq \\\n            -p 8 --validateMappings -o $1/${sample_name}_quant\ndone')


# ### Consolidate sample quantification to gene expression dataframe

# In[10]:


# PAO1
# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths_phage.PHAGE_QUANT

expression_phage_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_phage_df.head()


# In[11]:


# Save gene expression data
expression_phage_df.to_csv(paths_phage.PHAGE_GE, sep='\t')

