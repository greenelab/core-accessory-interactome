
# coding: utf-8

# # Process phage reference genome
# 
# Phage genomes were downloaded from NCBI GenBank using search keywords: [phage] AND [pseudomonas].
# 
# This search returned 1,950 samples (as of 15 December 2020)
# 
# By manual inspection, this download includes phages for other bacteria such as Samonella and E. Coli. This notebook removes those FASTA entries and saves the cleaned version.

# In[1]:


from Bio import SeqIO
from core_acc_modules import paths


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
SeqIO.write(cleaned_records, paths.PHAGE_REF, "fasta")

