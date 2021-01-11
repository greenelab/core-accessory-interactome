
# coding: utf-8

# # Build reference transcriptome
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

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

from core_acc_modules import paths


# In[2]:


# Get PAO1 index
get_ipython().system(' salmon index -t $paths.PAO1_REF -i $paths.PAO1_INDEX')


# In[3]:


# Get PA14 index
get_ipython().system(' salmon index -t $paths.PA14_REF -i $paths.PA14_INDEX')


# In[4]:


# Get phage index
get_ipython().system(' salmon index -t $paths.PHAGE_REF -i $paths.PHAGE_INDEX')

