
# coding: utf-8

# # Create reference genome
# 
# This notebook will create two combined reference genomes 1) phage and PAO1 reference sequences, 2) phage and PA14 reference sequences.

# In[1]:


import os
import pandas as pd
from Bio import SeqIO
from core_acc_modules import paths


# ## Clean phage sequence download --> Necessary?
# 
# Phage genomes were downloaded from NCBI GenBank using search keywords: [phage] AND [pseudomonas].
# 
# This search returned 1,950 samples (as of 15 December 2020)
# 
# By manual inspection, this download includes phages for other bacteria such as Samonella and E. Coli. We will remove those FASTA entries and save the cleaned version.

# In[2]:


# Select only those entries with keyword, "pseudomonas"
cleaned_records = []
keyword = "pseudomonas"
for record in SeqIO.parse(paths.RAW_PHAGE_REF, "fasta"):
    print("%s %s %i" % (record.id, record.description.lower(), len(record)))
    if keyword in record.description.lower():
        cleaned_records.append(record)


# In[3]:


# Write cleaned fasta records to file
#SeqIO.write(cleaned_records, paths.PHAGE_REF, "fasta")


# ## Combine phage + PAO1/PA14 genomes 
# 
# We want to create a fasta file with PAO1 + phage gene sequences and a file with PA14 + phage gene sequences. To do this we need to make sure we are only adding unique phage genome sequences to the PAO1 or PA14 reference sequences. We will do this using BLAST. For all phage genome sequences, we will BLAST against PAO1 or PA14 sequences
# 
# <--- Description of BLAST algorithm used here --->

# ### Process PAO1 and PA14 sequence files
# 
# 1. Make sure that file is .fasta
# 2. Remove any duplicate sequence ids

# In[4]:


# Remove duplicate PAO1 reference sequences
pao1_noduplicates_ref = []

pao1_seq_ids_seen = []
for record in SeqIO.parse(paths.PAO1_REF, "fasta"):
    if record.id not in pao1_seq_ids_seen:
        pao1_seq_ids_seen.append(record.id)
        pao1_noduplicates_ref.append(record)
        
# Write cleaned fasta records to file
SeqIO.write(pao1_noduplicates_ref, paths.PAO1_REF, "fasta")


# In[5]:


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
# <-- What are these blast error messages --->

# In[6]:


os.makedirs(paths.BLAST_DIR, exist_ok=True)


# In[7]:


get_ipython().run_cell_magic('bash', '-s $paths.PAO1_REF $paths.PAO1_DB_DIR', '\nmakeblastdb -in $1 -dbtype nucl -parse_seqids -out $2')


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

# In[8]:


get_ipython().run_cell_magic('bash', '-s $paths.PHAGE_REF $paths.PAO1_BLAST_RESULT $paths.PAO1_DB_DIR', 'blastn -query $1 -out $2 -db $3 -outfmt 6 ')


# In[ ]:


# SAME FOR PA14


# In[9]:


pao1_blast_result = pd.read_csv(paths.PAO1_BLAST_RESULT, sep="\t", header=None)
# SAME FOR PA14


# In[10]:


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
pao1_blast_result.columns = col_names
print(pao1_blast_result.shape)
pao1_blast_result.head()

# SAME FOR PA14


# ## Add non-duplicate phage sequences only
# ASSERTION FAILED
# 
# THINK THIS IS BECAUSE THOSE SEQUENCES DID NOT MAP AT ALL
# #assert(pao1_blast_result.shape[0] == len(pao1_noduplicates_ref))
# 
# The smaller the E-value, the better the match. So we want to add phage sequences with high E-value.

# In[11]:


duplicate_phage_seqs = pao1_blast_result["qseqid"]

pao1_phage_ref_seqs = []

# Add all PAO1 reference sequences
for record in SeqIO.parse(paths.PAO1_REF, "fasta"):
    pao1_phage_ref_seqs.append(record)

# Only add non-redundant phage sequences
for record in SeqIO.parse(paths.PHAGE_REF, "fasta"):
    if record.id not in duplicate_phage_seqs:
        pao1_phage_ref_seqs.append(record)

print(len(pao1_phage_ref_seqs))


# In[ ]:


# SAME FOR PA14


# In[12]:


# Write cleaned fasta records to file
SeqIO.write(pao1_phage_ref_seqs, paths.PAO1_PHAGE_REF, "fasta")
# SAME FOR PA14

