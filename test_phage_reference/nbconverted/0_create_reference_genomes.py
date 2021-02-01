#!/usr/bin/env python
# coding: utf-8

# # Create reference genome
# 
# This notebook will create two combined reference genomes 1) phage and PAO1 reference sequences, 2) phage and PA14 reference sequences.

# In[1]:


import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from core_acc_modules import paths


# ## Download phage genomes
# Phage genomes were downloaded from NCBI GenBank using search keywords: [phage] AND [pseudomonas].
# 
# This search returned 1,950 samples (as of 15 December 2020)

# ## Combine phage + PAO1/PA14 genomes 
# 
# We want to create a fasta file with PAO1 + phage gene sequences and a file with PA14 + phage gene sequences. To do this we need to make sure we are only adding unique phage genome sequences to the PAO1 or PA14 reference sequences. We will do this using BLAST. For all phage genome sequences, we will BLAST against PAO1 or PA14 sequences
# 
# **BLAST (basic local alignment search tool)** which performs comparisons between pairs of sequences, searching for regions of local similarity. 
# 
# 1. An initial search is done for a word of length "W" (W-mer) in the query sequence. 
# 2. Word hits are then extended in either direction in an attempt to generate a maximal segment pair (MSP)
# 3. The MSP score is computed based on the number of mismatches/matches, gaps
# 
# https://www.ncbi.nlm.nih.gov/books/NBK62051/

# ### Process PAO1 and PA14 sequence files
# 
# 1. Make sure that file is .fasta
# 2. Remove any duplicate sequence ids

# In[2]:


# Remove duplicate PAO1 reference sequences
pao1_noduplicates_ref = []

pao1_seq_ids_seen = []
for record in SeqIO.parse(paths.PAO1_REF, "fasta"):
    if record.id not in pao1_seq_ids_seen:
        pao1_seq_ids_seen.append(record.id)
        pao1_noduplicates_ref.append(record)
        
# Write cleaned fasta records to file
SeqIO.write(pao1_noduplicates_ref, paths.PAO1_REF, "fasta")


# In[3]:


# Remove duplicate PAO1 reference sequences
pa14_noduplicates_ref = []

pa14_seq_ids_seen = []
for record in SeqIO.parse(paths.PA14_REF, "fasta"):
    if record.id not in pa14_seq_ids_seen:
        pa14_seq_ids_seen.append(record.id)
        pa14_noduplicates_ref.append(record)
        
# Write cleaned fasta records to file
SeqIO.write(pa14_noduplicates_ref, paths.PA14_REF, "fasta")


# ### Create BLAST DB
# 
# Note: Got the following error message:
# 
# ```
# Error: (803.7) [makeblastdb] Blast-def-line-set.E.title
# Bad char [0xC2] in string at byte 52
# cds chromosome:3013928-3015481(+) name=fhpR ;locus_tag=PA2665;replicon_accession=NC_002516.2;product=Transcriptional activator of P. aeruginosa flavohemoglobin, FhpR
# ```
# 
# The error seems to be associated with non-acii characters present. https://github.com/tseemann/abricate/issues/77
# 
# **Fix:** Manually changed `name=fhpR ;` to `name=fhpR;`

# In[4]:


os.makedirs(paths.BLAST_DIR, exist_ok=True)


# In[5]:


get_ipython().run_cell_magic('bash', '-s $paths.PAO1_REF $paths.PAO1_DB_DIR', '\nmakeblastdb -in $1 -dbtype nucl -parse_seqids -out $2')


# In[6]:


get_ipython().run_cell_magic('bash', '-s $paths.PA14_REF $paths.PA14_DB_DIR', '\nmakeblastdb -in $1 -dbtype nucl -parse_seqids -out $2')


# ### BLAST phage sequences against PAO1/PA14 DB
# 
# http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# 
# | Column header | Description |
# | --- | --- | 
# | qseqid | query (e.g., unknown gene) sequence id |
# | sseqid | subject (e.g., reference genome) sequence id |
# | pident | percentage of identical matches |
# | length | alignment length (sequence overlap) |
# | mismatch | number of mismatches |
# | gapopen | number of gap openings |
# | qstart | start of alignment in query |
# | qend | end of alignment in query |
# | sstart | start of alignment in subject
# | send | end of alignment in subject |
# | evalue | expect value |
# | bitscore | bit score|
# 
# E-value: expected number of chance alignments; the smaller the E-value, the better the match.

# In[7]:


get_ipython().run_cell_magic('bash', '-s $paths.PHAGE_REF $paths.PAO1_BLAST_RESULT $paths.PAO1_DB_DIR', 'blastn -query $1 -out $2 -db $3 -outfmt 6 ')


# In[8]:


get_ipython().run_cell_magic('bash', '-s $paths.PHAGE_REF $paths.PA14_BLAST_RESULT $paths.PA14_DB_DIR', 'blastn -query $1 -out $2 -db $3 -outfmt 6 ')


# In[9]:


pao1_blast_result = pd.read_csv(paths.PAO1_BLAST_RESULT, sep="\t", header=None)
pa14_blast_result = pd.read_csv(paths.PA14_BLAST_RESULT, sep="\t", header=None)


# In[10]:


# Add column names described above
col_names = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
    
]


# In[11]:


# BLAST results for PAO1
pao1_blast_result.columns = col_names
print(pao1_blast_result.shape)
print(pao1_blast_result["evalue"].max())
pao1_blast_result.head()


# In[12]:


# BLAST results for PA14
pa14_blast_result.columns = col_names
print(pa14_blast_result.shape)
print(pa14_blast_result["evalue"].max())
pa14_blast_result.head()


# ## Add non-duplicate phage sequences only
# The smaller the E-value, the better the match. So we want to add phage sequences with high E-value (i.e. evalue > 0.05) or phage sequences that are not on this table at all, because there was no hit found.
# 
# In this case, there were no sequences with high evalues so we are treating all the sequences in the above table as duplicates

# In[13]:


duplicate_pao1_phage_seq_ids = np.unique(pao1_blast_result["qseqid"].values)

pao1_phage_ref_seqs = []

# Add all PAO1 reference sequences
for record in SeqIO.parse(paths.PAO1_REF, "fasta"):
    pao1_phage_ref_seqs.append(record)

num_pao1_seqs = len(pao1_phage_ref_seqs)
print(num_pao1_seqs)

# Only add non-redundant phage sequences
num_phage_seqs = 0
for record in SeqIO.parse(paths.PHAGE_REF, "fasta"):
    num_phage_seqs += 1
    if record.id not in duplicate_pao1_phage_seq_ids:
        pao1_phage_ref_seqs.append(record)

print(num_phage_seqs)
print(len(duplicate_pao1_phage_seq_ids))
print(len(pao1_phage_ref_seqs))

# Check length
assert len(pao1_phage_ref_seqs) == (num_pao1_seqs + num_phage_seqs - len(duplicate_pao1_phage_seq_ids))


# In[14]:


duplicate_pa14_phage_seq_ids = np.unique(pa14_blast_result["qseqid"].values)

pa14_phage_ref_seqs = []

# Add all PAO1 reference sequences
for record in SeqIO.parse(paths.PA14_REF, "fasta"):
    pa14_phage_ref_seqs.append(record)

num_pa14_seqs = len(pa14_phage_ref_seqs)
print(num_pa14_seqs)

# Only add non-redundant phage sequences
for record in SeqIO.parse(paths.PHAGE_REF, "fasta"):
    if record.id not in duplicate_pa14_phage_seq_ids:
        pa14_phage_ref_seqs.append(record)

print(num_phage_seqs)
print(len(duplicate_pa14_phage_seq_ids))
print(len(pa14_phage_ref_seqs))

# Check length
assert len(pa14_phage_ref_seqs) == (num_pa14_seqs + num_phage_seqs - len(duplicate_pa14_phage_seq_ids))


# In[15]:


# Write cleaned fasta records to file
SeqIO.write(pao1_phage_ref_seqs, paths.PAO1_PHAGE_REF, "fasta")
SeqIO.write(pa14_phage_ref_seqs, paths.PA14_PHAGE_REF, "fasta")

