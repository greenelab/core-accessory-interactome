
# coding: utf-8

# # Differential expression validation
# This notebook performs a differential expression (DE) analysis comparing PAO1 samples vs PA14 samples. We can compare our results with those published in the literature as an additional step to validate that our RNA-seq processing are reasonable.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from core_acc_modules import utils, paths
from rpy2.robjects import pandas2ri
pandas2ri.activate()


# ### Download data for validation
# Data from [Sana et. al](https://jb.asm.org/content/201/21/e00362-19) found ~ 2K DEGs between 2 strains where QS genes were DEGs.

# In[2]:


# Download sra data files
get_ipython().system(' prefetch --option-file $paths.SRA_ACC_TEST')


# In[3]:


get_ipython().run_cell_magic('bash', '', 'mkdir $paths.FASTQ_TEST_DIR\nfastq-dump $paths.SRA_DIR/SRR8486287.sra --outdir $paths.FASTQ_TEST_DIR/\nfastq-dump $paths.SRA_DIR/SRR8486288.sra --outdir $paths.FASTQ_TEST_DIR/\nfastq-dump $paths.SRA_DIR/SRR8486289.sra --outdir $paths.FASTQ_TEST_DIR/\nfastq-dump $paths.SRA_DIR/SRR8486290.sra --outdir $paths.FASTQ_TEST_DIR/')


# In[4]:


get_ipython().run_cell_magic('bash', '-s $paths.PAO1_QUANT_TEST $paths.FASTQ_TEST_DIR $paths.PAO1_INDEX', 'mkdir $1\n\nfor FILE_PATH in $2/*;\ndo\n\n# get file name\nsample_name=`basename ${FILE_PATH}`\n\n# remove extension from file name\nsample_name="${sample_name%_*}"\n\n# get base path\nbase_name=${FILE_PATH%/*}\n\necho "Processing sample ${sample_name}"\n\nsalmon quant -i $3 -l A \\\n            -r ${base_name}/${sample_name} \\\n            -p 8 --validateMappings -o $1/${sample_name}_quant\ndone')


# In[5]:


# Get raw read counts using PAO1 reference
# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PAO1_QUANT_TEST

expression_data = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["NumReads"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

# Map gene ids to gene names
pao1_fasta_file = paths.PAO1_REF

seq_id_to_gene_id_pao1 = utils.dict_gene_num_to_ids(pao1_fasta_file)

expression_data.rename(mapper=seq_id_to_gene_id_pao1, axis="columns", inplace=True)

expression_data.head()


# In[6]:


new_index = [name.split(".")[0] for name in list(expression_data.index)]
expression_data.index = new_index
expression_data.head()


# ### Process data
# 1. Get core genes
# 2. Round read counts to integer value

# In[7]:


# Get mapping between PAO1 and PA14 genes using PAO1 reference
gene_annot_file = paths.GENE_PAO1_ANNOT
gene_mapping_pao1 = utils.get_pao1_pa14_gene_map(gene_annot_file, 'pao1')
gene_annot_file = paths.GENE_PA14_ANNOT
gene_mapping_pa14 = utils.get_pao1_pa14_gene_map(gene_annot_file, 'pa14')

core_pao1_genes, core_pa14_genes = utils.get_core_genes(gene_mapping_pao1,
                                                        gene_mapping_pa14,
                                                        False)
print(f"Number of PAO1 core genes: {len(core_pao1_genes)}")
print(f"Number of PA14 core genes: {len(core_pa14_genes)}")

core_pao1_genes = set(core_pao1_genes) - set(["PA4215", "PA4214","PA4213"])

expression_data = expression_data.reindex(columns=core_pao1_genes)
print(expression_data.shape)
expression_data.head()


# In[8]:


# Convert values to integers for DE analysis
# Not sure why the "Numreads" values are floats
expression_data = expression_data.astype(int)


# In[9]:


# Save file
expression_data.to_csv(paths.PAO1_GE_DE, sep='\t')


# ### Examine gene expression
# 1. Look at consistency of gene expression within PAO1 samples, within PA14 samples
# 2. What does average PAO1 vs averge PA14 look like? To give us an expectation for our DESeq analysis.

# In[10]:


expression_data = pd.read_csv(paths.PAO1_GE_DE, sep="\t", index_col=0)
expression_data.head()


# In[11]:


# Group samples as PAO1 or PA14 based on experiment metadata
sample_annot_file = paths.SAMPLE_ANNOT_TEST
pao1_ids, pa14_ids = utils.get_sample_grps(sample_annot_file)


# In[12]:


# Split expression by genotype
pao1_expression = expression_data.loc[pao1_ids]
pa14_expression = expression_data.loc[pa14_ids]


# In[13]:


# Get within sample correlation
pao1_corr = pao1_expression.T.corr()
pa14_corr = pa14_expression.T.corr()


# In[14]:


ax = sns.heatmap(
    pao1_corr, 
    vmin=-1, vmax=1, center=0,
    cmap=sns.diverging_palette(20, 220, n=200),
    square=True
)
ax.set_title("PAO1 sample correlation")


# In[15]:


ax = sns.heatmap(
    pa14_corr, 
    vmin=-1, vmax=1, center=0,
    cmap=sns.diverging_palette(20, 220, n=200),
    square=True
)
ax.set_title("PA14 sample correlation")


# In[16]:


# Get mean expression for pao1 and pa14
mean_pao1 = pao1_expression.mean()
mean_pa14 = pa14_expression.mean()

pao1_v_pa14_df = pd.DataFrame(data={'pao1_mean': mean_pao1.values,
                                   'pa14_mean': mean_pa14.values},
                             index=pao1_expression.columns)

pao1_v_pa14_df.head()


# In[17]:


sns.scatterplot(data=pao1_v_pa14_df, x='pao1_mean', y='pa14_mean')
plt.plot([0,140000],[0,140000])


# In[18]:


# Rough calculation of the number of genes that differ
# This is roughly what we would expect to get from the DE analysis

std_x_eq_y = np.std(abs(pao1_v_pa14_df['pao1_mean']-pao1_v_pa14_df['pa14_mean']))
gene_differences = pao1_v_pa14_df[abs(pao1_v_pa14_df['pao1_mean']-pao1_v_pa14_df['pa14_mean']) > std_x_eq_y]
print(gene_differences.shape)
genes_found_from_GE = list(gene_differences.index)
gene_differences.head()


# **Observations:**
# * Looks like there is consistent gene expression patterns within sample-type (i.e. Both PAO1 samples have a similar gene expression profile, similarly for both PA14 samples) as expected
# * Comparing the mean expression of PAO1 and PA14 we see that there are ~200 genes changed. This gives us some indication about what to expect for our DESeq analysis. However, we shouldn't expect the numbers to align because we are using different methods -- above we are comparing the raw gene expression values and looking for a threshold difference; below we DESeq fits the negative binomial model to the data and performs hypothesis testing to determine if there is a difference between the groups of samples.  

# ### Differential expression analysis

# In[19]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Run one time\n#if (!requireNamespace("BiocManager", quietly = TRUE))\n#    install.packages("BiocManager")\n#BiocManager::install("DESeq2")')


# In[20]:


get_ipython().run_cell_magic('R', '', '# Load the DESeq2 library\nlibrary("DESeq2")')


# In[21]:


# Files to load into DE analysis (R)
metadata_file = str(paths.SAMPLE_ANNOT_TEST)
expression_data_file = str(paths.PAO1_GE_DE)
out_file = str(paths.DE_STATS)


# In[22]:


# Check ordering of sample ids
utils.check_sample_ordering(expression_data_file, metadata_file)


# In[23]:


get_ipython().run_cell_magic('R', '-i metadata_file -i expression_data_file -i out_file', "\nsource('../core_acc_modules/DE_analysis.R')\n\nget_DE_stats_DESeq(metadata_file,\n                   expression_data_file,\n                   out_file)")


# In[24]:


# Read in DE stats file
DE_stats = pd.read_csv(paths.DE_STATS, sep='\t', header=0, index_col=0)
print(DE_stats.shape)
DE_stats.head()


# In[25]:


# Volcano plot
DE_stats["-log10(padj)"] = -np.log10(DE_stats["padj"])
sns.scatterplot(data=DE_stats, x="log2FoldChange", y="-log10(padj)")


# ### Compare our DE results with publication

# In[26]:


# Get number of DEGs
selected_DE_stats = DE_stats[(DE_stats['padj']<0.01)]

print(selected_DE_stats.shape)
selected_DE_stats.head()


# In[27]:


# Compare our findings against Sana et. al.
degs_sana = ["PA3431", "PA3432", "PA1244", "PA4685"]
for gene in degs_sana:
    if gene in list(selected_DE_stats.index):
        print(gene)


# In[28]:


# Compare genes whose mean gene expression differed between PAO1 and PA14
# with those genes found using DESeq
top_degs_by_padj = selected_DE_stats.index
len(set(top_degs_by_padj).intersection(genes_found_from_GE))/len(genes_found_from_GE)


# **Conclusions:**
# 
# Our DE analysis found ~1.1K significantly differentially expressed genes
# 
# (Check 1): The DEGs identified using DESeq (\~1.1K genes) is fairly consistent with the genes that were 1 standard deviation outside the correlation threshold (\~200 genes) -- there was a 75% overlap in these gene sets. This very roughly validates that DESeq is working as expected. I wouldn't expect the numbers to be this different though.
# 
# (Check 2) The number of DEGs identified (\~1.1K genes) using DESeq is fairly consistent with the number of differentially expressed genes found in [Sana et. al](https://jb.asm.org/content/201/21/e00362-19) (\~2K genes). We also spot checked specific genes that were found. We found the 4 genes highlighted in the Sana et. al. publication, including the main qsIA gene (PA1244) that the paper found to be more highly expressed in PAO1 vs PA14. Difference are likely due to differences in the package used. 
# 
# Approach used in [Sana et. al](https://jb.asm.org/content/201/21/e00362-19) found ~ 2K DEGs between 2 strains where QS genes were DEGs:
# ```
# Illumina reads were mapped to the P. aeruginosa genome PAO1 (GenBank accession number AE004091.2 [61]) and PA14 (GenBank accession number NC_008463.1 [33]) by Bowtie (version Bowtie1 v0.12.9 [62]). Data were normalized by reads per kilobase per million (RPKM) and filtered to the 5,263 orthologous genes conserved between P. aeruginosa strains PA14 and PAO1. Two biological replicates were performed per condition. Differential
# expression analysis was analyzed using the Bioconductor package NOISeq version 2.22.1 (64), a nonparametric approach suitable for lowly replicated data, and using a q value of 0.99 for strong control of
# false positives
# ```
