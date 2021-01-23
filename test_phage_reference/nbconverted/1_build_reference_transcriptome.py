
# coding: utf-8

# # Build reference transcriptome
# 
# The first step to using Salmon to get gene expression data is to build a reference transcriptome to align our sample reads against.
# 
# **Goal:** To create an index to evaluate sequences for all possible unique sequences of length k in the target transcriptome
# 
# **Input:** 
# 
# Target transcriptome. This transcriptome is given in the form of a multi-FASTA file, with each entry providing the sequence of a transcript. For prokaryotes, transcripts and genes have a more 1-1 mapping so we're using genes for our reference transcriptome and so we don't need to use tximport to map transcript quants to genes.
# 
# *What is the relationship between a gene vs transcript?* A gene is a bounded region of a chromosome within which transcription occurs. The gene (DNA sequence) is transcribed into RNA transcript. In bacteria, these RNA transcripts act as mRNA that can be translated into protein. In eukaryotes, the transcript RNA is pre-mRNA and must undergo additional processing (post-transcriptional modifications) before it can be translated. This processing includes addition of a protective cap and tail, splicing to remove introns. So genes can have multiple mRNA (which can encode different proteins) through the process of alternative splicing, where fragments of the pre-mRNA are assembled in different ways.
# 
# **Output:**
# 
# Quasi-index over the reference transcriptome, which is a structure that salmon uses to quasi-map RNA-seq reads during quantification. The index is a signature for each transcript in the reference transcriptome
# 
# The index contains:
# * Suffix array (SA) of the reference transcriptome. There is a suffix array per transcript in the reference, containing a sorted array of all the suffixes of each transcript
# * A hash table mapping each k-mer occurring in the reference transcriptome to its location in SA
# * https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/08_salmon.html
# 
# **Command:**
# 
# `> ./bin/salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31`
# * k = minimum acceptable length for valid matches. So a smaller k might improve sensitivity. They found that k=31 works well for reads of 75bp or longer. 
# * decoys = Set of target transcript ids that will not appear in the quantification. This decoy transcriptome is meant to mitigate potential spurious mapping of reads that actually arise from some unannotated genomic locus that is sequence-similar to an annotated transcriptome.
# 
# 
# **Note:** Here we are using [Salmon](https://combine-lab.github.io/salmon/) version 0.11.2 to be consistent with the version running on Dartmouth's computing cluster (Discovery), which is where the data will be processed. 

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

from core_acc_modules import paths


# In[2]:


# Get PAO1 index
get_ipython().system(' salmon index -t $paths.PAO1_PHAGE_REF -i $paths.PAO1_PHAGE_INDEX')


# In[3]:


# Get PA14 index
get_ipython().system(' salmon index -t $paths.PA14_PHAGE_REF -i $paths.PA14_PHAGE_INDEX')


# **Thoughts based on output:**
# * How does this handle full genomes for phages? Each entry will be much longer in length than Salmon expects so a warning message is output: `Entry with header [NC_028999.1] was longer than 200000 nucleotides.  This is probably a chromosome instead of a transcript.` Is Salmon including these entries?
# * When building the index I'm getting that PAO1 34 duplicates are removed, PA14 37 duplicates are removed, Phage  391 duplicates removed. Are these duplicates expected?
