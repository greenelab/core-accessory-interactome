
# coding: utf-8

# # Download and process SRA data

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
from pathlib import Path
import pandas as pd
import numpy as np
import umap
import seaborn as sns
import matplotlib.pyplot as plt
from core_acc_modules import utils, paths
from plotnine import (ggplot,
                      labs,  
                      geom_point,
                      aes, 
                      ggsave,
                      theme,
                      theme_bw,
                      scale_color_manual,
                      guides, 
                      guide_legend,
                      element_blank,
                      element_text,
                      element_rect,
                      element_line,
                      coords)

np.random.seed(123)


# ### Setup SRA toolkit

# In[2]:


# Download latest version of compiled binaries of NCBI SRA toolkit 
if not os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
    get_ipython().system(' wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"')


# In[3]:


# Extract tar.gz file 
if os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
    get_ipython().system(' tar -xzf sratoolkit.current-centos_linux64.tar.gz')

# add binaries to path using export path or editing ~/.bashrc file
get_ipython().system(' export PATH=$PATH:sratoolkit.2.10.7-centos_linux64/bin')

# Now SRA binaries added to path and ready to use


# ### Download SRA data
# 
# Two SRA projects were selected:
# * [PRJNA633671](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA633671) 13 PAO1 samples isolated from pig burn wound after 3 (3), 14(5), 28(5) days post infection
# * [PRJNA491911](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA491911) 15 PA14 samples isolated from pig burn wound after 3 (5), 14 (5), 28 (5) days post infection

# In[4]:


# Download sra data files
get_ipython().system(' prefetch --option-file data/metadata/sra_acc.txt ')


# ### Get FASTQ files associated with SRA downloads
# 
# The fastq files store the RNA-seq results, including: sequencing and quality scores for each base call.
# 
# Here is a nice blog to explain how to read fastq files: https://thesequencingcenter.com/knowledge-base/fastq-files/
# 
# The fastq files gives the sequence of a read at a given location. Our goal is to map these reads to a reference genome so that we can quantify the number of reads that are at a given location, to determine the level of expression.

# In[5]:


get_ipython().system('mkdir $paths.FASTQ_DIR')


# In[6]:


get_ipython().system('fastq-dump $paths.SRA_DIR/* --split-files --outdir $paths.FASTQ_DIR/')


# ### Obtain a transcriptome and build an index
# 
# Here we are using [Salmon](https://combine-lab.github.io/salmon/)
# 
# **Input:**
# * Target transcriptome
# * This transcriptome is given to Salmon in the form of a (possibly compressed) multi-FASTA file, with each entry providing the sequence of a transcript
# * DNA sequences (genes) that get transcribed as mRNA (transcripts). Should we look at reference transcipts or genes?
# * We downloaded the `GENE DNA` file for `Pseudomonas aeruginosa PAO1 (Reference)` and `Pseudomonas aeruginosa UCBPP-PA14` from  http://www.pseudomonas.com/strain/download
# 
# **!!!** Georgia is using PAO1 reference from ensembl: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-47/fasta/bacteria_13_collection/
# Don't see an equivalent for PA14.
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


# Get PAO1 index
get_ipython().system(' salmon index -t $paths.PAO1_REF -i $paths.PAO1_INDEX')


# In[8]:


# Get PA14 index
get_ipython().system(' salmon index -t $paths.PA14_REF -i $paths.PA14_INDEX')


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


get_ipython().run_cell_magic('bash', '-s $paths.PAO1_QUANT $paths.FASTQ_DIR $paths.PAO1_INDEX', 'mkdir $1\n\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${sample_name}"\n\nsalmon quant -i $3 -l A \\\n            -1 ${base_name}/${sample_name}_1.fastq \\\n            -2 ${base_name}/${sample_name}_2.fastq \\\n            -p 8 --validateMappings -o $1/${sample_name}_quant\ndone')


# #### Get quants using PA14 reference

# In[ ]:


get_ipython().run_cell_magic('bash', '-s $paths.PA14_QUANT $paths.FASTQ_DIR $paths.PA14_INDEX', 'mkdir $1\n\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${sample_name}"\n\nsalmon quant -i $3 -l A \\\n            -1 ${base_name}/${sample_name}_1.fastq \\\n            -2 ${base_name}/${sample_name}_2.fastq \\\n            -p 8 --validateMappings -o $1/${sample_name}_quant\ndone')


# ### Consolidate sample quantification to gene expression dataframe

# In[2]:


# PAO1
# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PAO1_QUANT

expression_pao1_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_pao1_df.head()


# In[3]:


# PA14
data_dir = paths.PA14_QUANT

expression_pa14_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

expression_pa14_df.head()


# In[4]:


# Map gene ids to gene names
pa14_fasta_file = paths.PA14_REF
pao1_fasta_file = paths.PAO1_REF

seq_id_to_gene_id_pao1 = utils.dict_gene_num_to_ids(pao1_fasta_file)
seq_id_to_gene_id_pa14 = utils.dict_gene_num_to_ids(pa14_fasta_file)

expression_pao1_df.rename(mapper=seq_id_to_gene_id_pao1, axis="columns", inplace=True)
expression_pa14_df.rename(mapper=seq_id_to_gene_id_pa14, axis="columns", inplace=True)


# In[5]:


# Save gene expression data
expression_pao1_df.to_csv(paths.PAO1_GE, sep='\t')
expression_pa14_df.to_csv(paths.PA14_GE, sep='\t')


# ### Quick validation
# Here we want to validate that we've processed the samples correctly using Salmon.

# In[6]:


# Get PAO1 core and accessory genes
gene_mapping = utils.get_pao1_pa14_gene_map(paths.GENE_PAO1_ANNOT, 'pao1')
core_genes = utils.get_core_genes(gene_mapping)
acc_genes = list(set(expression_pao1_df.columns) - set(core_genes))

print(len(core_genes))
print(len(acc_genes))


# In[7]:


assert("PA0053" in acc_genes)


# In[8]:


# Load in sample annotation file
sample_annot_file = paths.SAMPLE_ANNOT
pao1_ids, pa14_ids = utils.get_sample_grps(sample_annot_file)


# In[9]:


# Examine PA14 samples in PAO1-specific genes (PAO1 reference)
pa14_samples_pao1_genes_pao1_ref = expression_pao1_df.loc[pa14_ids,acc_genes]
pa14_samples_pao1_genes_pao1_ref_mean = pa14_samples_pao1_genes_pao1_ref.mean()
pa14_samples_pao1_genes_pao1_ref_mean.isna().any()


# In[10]:


# Examine PA14 samples in core genes
pa14_samples_core_genes_pao1_ref = expression_pao1_df.loc[pa14_ids,core_genes]
pa14_samples_core_genes_pao1_ref_mean = pa14_samples_core_genes_pao1_ref.mean()
pa14_samples_core_genes_pao1_ref_mean.isna().any()
pa14_samples_core_genes_pao1_ref_mean[pa14_samples_core_genes_pao1_ref_mean.isna()]


# In[11]:


# Plot
sns.set()

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(6,6))

# Distribution plot for core genes
sns.distplot(pa14_samples_core_genes_pao1_ref_mean.values,
             label='PA14 samples core genes',
             color='red',
             #bins=bins_expression,
             kde=False,
             ax=axes[0]
            )

sns.distplot(pa14_samples_pao1_genes_pao1_ref_mean.values, 
             label='PA14 samples PAO1 specific genes', 
             color='blue',
             #bins=bins_expression,
             kde=False,
             ax=axes[1]
            )

fig.xlim=(0,10)
plt.suptitle('Histogram of mean gene expression of PA14 samples (PAO1 reference)',
            fontsize=16)
fig.text(0.5, 0.01, 'Mean gene expression', ha='center', fontsize=14)
fig.text(0.01, 0.5, 'Count', ha='center', rotation=90, fontsize=14)
plt.tight_layout(pad=0.4, 
                 w_pad=0.5,
                 h_pad=1.0,
                 rect=[0, 0.03, 1, 0.95])


# **Takeaway**:
# The plot above is taking all PA14 samples and looking at the distribution of mean gene expression(across samples) in two cases: (blue) mean gene expression for PAO1-specific genes (i.e. genes absent in PA14 strains) and (red) mean gene expression for core genes (i.e. genes shared by both PAO1 and PA14 strains). 
# 
# If we processed the data correctly, we'd expect that the mean expression of PAO1-specific genes (shown in blue) (ie. those PAO1 genes that do not have a PA14 homolog) have 0 expression in PA14 samples. We see that most PAO1-specific genes do have 0 expression. In comparison the mean expression of the core genes are mainly nonzero. 

# #### Visualize clustering of gene expression

# In[12]:


# Embed expression data into low dimensional space
model = umap.UMAP(random_state=123).fit(expression_pao1_df)
pao1_encoded = model.transform(expression_pao1_df)

pao1_encoded_df = pd.DataFrame(data=pao1_encoded,
                               index=expression_pao1_df.index,
                               columns=['1','2'])

# Add label
pao1_encoded_df['genotype'] = 'PAO1'
pao1_encoded_df.loc[pa14_ids,'genotype'] = 'PA14'

pao1_encoded_df.head()


# In[13]:


# Embed expression data into low dimensional space
model = umap.UMAP(random_state=123).fit(expression_pa14_df)
pa14_encoded = model.transform(expression_pa14_df)

pa14_encoded_df = pd.DataFrame(data=pa14_encoded,
                               index=expression_pa14_df.index,
                               columns=['1','2'])

# Add label
pa14_encoded_df['genotype'] = 'PAO1'
pa14_encoded_df.loc[pa14_ids,'genotype'] = 'PA14'

pa14_encoded_df.head()


# In[14]:


# Plot PAO1
fig = ggplot(pao1_encoded_df, aes(x='1', y='2'))
fig += geom_point(aes(color='genotype'), alpha=0.5)
fig += labs(x ='UMAP 1',
            y = 'UMAP 2',
            title = 'RNA-seq expression using PAO1 reference')
fig += theme_bw()
fig += theme(
    legend_title_align = "center",
    plot_background=element_rect(fill='white'),
    legend_key=element_rect(fill='white', colour='white'), 
    legend_title=element_text(family='sans-serif', size=15),
    legend_text=element_text(family='sans-serif', size=12),
    plot_title=element_text(family='sans-serif', size=15),
    axis_text=element_text(family='sans-serif', size=12),
    axis_title=element_text(family='sans-serif', size=15)
    )
fig += guides(colour=guide_legend(override_aes={'alpha': 1}))

print(fig)


# In[15]:


# Plot PA14
fig = ggplot(pa14_encoded_df, aes(x='1', y='2'))
fig += geom_point(aes(color='genotype'), alpha=0.5)
fig += labs(x ='UMAP 1',
            y = 'UMAP 2',
            title = 'RNA-seq expression using PA14 reference')
fig += theme_bw()
fig += theme(
    legend_title_align = "center",
    plot_background=element_rect(fill='white'),
    legend_key=element_rect(fill='white', colour='white'), 
    legend_title=element_text(family='sans-serif', size=15),
    legend_text=element_text(family='sans-serif', size=12),
    plot_title=element_text(family='sans-serif', size=15),
    axis_text=element_text(family='sans-serif', size=12),
    axis_title=element_text(family='sans-serif', size=15)
    )
fig += guides(colour=guide_legend(override_aes={'alpha': 1}))

print(fig)


# **Takeaway:**
# This plot is showing the clustering of samples using both the PAO1 reference transcriptome and the PA14 reference transcriptome. The plot shows that the samples clustering by genotype, as expected.
