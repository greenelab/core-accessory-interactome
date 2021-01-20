
# coding: utf-8

# # Align and quantify samples
# 
# The second step is to align our samples against our built index and quantify reads. This notebook aligns a pilot set of samples against the PAO1, PA14 and and PAO1 reference genomes. 
# 
# **Input:**
# * Index of reference transcriptome
# * FASTQ of experimental samples
# 
# **Output:**
# * For each query (fragment), the reference sequences (transcripts), strand and position from which the query may have likely originated. In many cases, this mapping information is sufficient for downstream analysis like transcript quantification
# * Each sample will have a quantification file (called quant.sf). The TSV file will contain the name (Name) of each transcript, its length (Length), effective length (EffectiveLength), and its abundance in terms of Transcripts Per Million (TPM) and estimated number of reads (NumReads) originating from this transcript.
#   * The first two columns are self-explanatory, the name of the transcript and the length of the transcript in base pairs (bp).
#   * The effective length represents the various factors that effect the length of transcript (i.e degradation, technical limitations of the sequencing platform)
#   * Salmon outputs ‘pseudocounts’ which predict the relative abundance of different isoforms in the form of three possible metrics (KPKM, RPKM, and TPM). TPM (transcripts per million) is a commonly used normalization method as described in [1] and is computed based on the effective length of the transcript.
#   * Estimated number of reads (an estimate of the number of reads drawn from this transcript given the transcript’s relative abundance and length)
# 
# **Algorithm:**
# * Given the index and a set of sequenced reads, the quant command quasi-maps the reads and uses the resulting mapping information to estimate transcript abundances.
# * Quasi-map is a way to map sequenced fragments (single or paired-end reads) to a target transcriptome. Quasi-mapping produces what we refer to as fragment mapping information. In particular, it provides, for each query (fragment), the reference sequences (transcripts), strand and position from which the query may have likely originated. In many cases, this mapping information is sufficient for downstream analysis like quantification.
# 
# Basic steps:
# 1. Scan read until k-mer appears in the hash table
# 2. Look up all SA intervals (reference suffixes) containing that k-mer
# 3. Maximum mappable prefix (MMP) finds the longest read sequence that exactly matches the reference suffix 
# 4. Determine the next informative position (NIP) by a k-mer skipping approach
# 5. Repeat until the end of the read
# 6. Reports all MMPs that read intersected with
# 
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/08_salmon.html
# 
# **Command:**
# 
# Command and parameters:
# 
# `> ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant`
# 
# * libtype = Determines how the reads should be interpreted including the relative orientation of paired ends (inward, outward, matching) and strandedness (stranded=specify if read1 comes from forward or reverse strand, unstranded). Currently using “A” which lets Salmon automatically decide
#   * https://salmon.readthedocs.io/en/latest/salmon.html
#   * (pg 38) https://buildmedia.readthedocs.org/media/pdf/salmon/stable/salmon.pdf
# * What does strand bias in output mean? Strand_bias is such that a value of 0.5 means there is no bias (i.e. half of the fragments start with read 1 on the forward strand and half start with read 1 on the reverse complement strand). 
#   * Based on lib_format_counts.json file, the bias is very close to 0.5, just above
#   * https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/08_salmon.html
# 
# 

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
import shutil
import pandas as pd
import numpy as np
from core_acc_modules import paths

np.random.seed(123)


# ### Quantify gene expression
# Now that we have our index built and all of our data downloaded, we’re ready to quantify our samples
# 
# **For each sample we have read counts per gene (where the genes are based on the reference gene file provided above).** 
# 
# Note about TPM calculation:
# * For sample A, transcript X will have read count
# * Reads per kilobase (RPK) = read count/length of the transcript
# * Per million scaling factor = sum(read count/length across all samples for that transcript)/1M
# * TPM = RPK/per million scaling factor
# * TPM will depend on the scaling factor. If the number of mapped reads is very low then scale factor will be very low and so any trancript that mapped will be increased to create outliers since we’re dividing by a small scale factor
# 
# **How are paired-end vs single-end read samples treated?**
# * The multiple files are treated as a single library, meaning?? In this case we are not specifying 
# 
# How are we reading the fastq files to be quantified? Are we reading in all files if multipler are provided?

# #### Get quants using PAO1 reference

# In[2]:


os.makedirs(paths.PAO1_QUANT, exist_ok=True)


# In[3]:


get_ipython().run_cell_magic('bash', '-s $paths.PAO1_QUANT $paths.FASTQ_DIR $paths.PAO1_INDEX', '\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${base_name}/${sample_name}/*"\n\nsalmon quant -i $3 \\\n             -l A \\\n             -r ${base_name}/${sample_name}/* \\\n             -o $1/${sample_name}_quant\ndone')


# #### Get quants using PA14 reference

# In[4]:


os.makedirs(paths.PA14_QUANT, exist_ok=True)


# In[5]:


get_ipython().run_cell_magic('bash', '-s $paths.PA14_QUANT $paths.FASTQ_DIR $paths.PA14_INDEX', '\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${base_name}/${sample_name}/*"\n\nsalmon quant -i $3 \\\n             -l A \\\n             -r ${base_name}/${sample_name}/* \\\n             -o $1/${sample_name}_quant\ndone')


# #### Get quants using phage reference

# In[6]:


os.makedirs(paths.PHAGE_QUANT, exist_ok=True)


# In[7]:


get_ipython().run_cell_magic('bash', '-s $paths.PHAGE_QUANT $paths.FASTQ_DIR $paths.PHAGE_INDEX', '\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${base_name}/${sample_name}/*"\n\nsalmon quant -i $3 \\\n             -l A \\\n             -r ${base_name}/${sample_name}/* \\\n             -o $1/${sample_name}_quant\ndone')


# ### Consolidate sample quantification to gene expression dataframe

# In[ ]:


# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PAO1_QUANT

expression_pao1_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_pao1_df.head()


# In[ ]:


# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PA14_QUANT

expression_pa14_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_pa14_df.head()


# In[ ]:


# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PHAGE_QUANT

expression_phage_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_phage_df.head()


# In[ ]:


# Save gene expression data
expression_pao1_df.to_csv(paths.PAO1_GE, sep='\t')
expression_pa14_df.to_csv(paths.PA14_GE, sep='\t')
expression_phage_df.to_csv(paths.PHAGE_GE, sep='\t')

