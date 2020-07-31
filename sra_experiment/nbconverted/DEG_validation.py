
# coding: utf-8

# # Differential expression validation
# This notebook performs a differential expression (DE) analysis comparing PAO1 samples vs PA14 samples. We can compare our results with those published in the literature as an additional step to validate that our RNA-seq processing are reasonable.

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd
import seaborn as sns
import numpy as np
from core_acc_modules import utils, paths
from rpy2.robjects import pandas2ri
pandas2ri.activate()


# In[2]:


"""# Get raw read counts using PAO1 reference
# Read through all sample subdirectories in quant/
# Within each sample subdirectory, get quant.sf file
data_dir = paths.PAO1_QUANT

expression_data = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["NumReads"].
    rename(file.parent.name.split("_")[0]) 
    for file in data_dir.rglob("*/quant.sf"))    

# Map gene ids to gene names
pao1_fasta_file = paths.PAO1_REF

seq_id_to_gene_id_pao1 = utils.dict_gene_num_to_ids(pao1_fasta_file)

expression_data.rename(mapper=seq_id_to_gene_id_pao1, axis="columns", inplace=True)

expression_data.head()"""


# In[3]:


# Load gene expression using PAO1 reference
expression_data = pd.read_csv(paths.PAO1_GE, sep='\t', header=0, index_col=0)
print(expression_data.shape)
expression_data.head()


# ### Get core genes

# In[4]:


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


# In[5]:


# Convert values to integers for DE analysis
# Not sure why the "Numreads" values are floats
#expression_data = expression_data.astype(int)


# In[6]:


# Save file
expression_data.to_csv(paths.PAO1_GE_DE, sep='\t')


# ### Differential expression analysis

# In[7]:


get_ipython().run_cell_magic('R', '', '# Select 59\n# Run one time\nif (!requireNamespace("BiocManager", quietly = TRUE))\n    install.packages("BiocManager")\n#BiocManager::install("DESeq2")\n#BiocManager::install("limma")')


# In[8]:


get_ipython().run_cell_magic('R', '', '# Load the DESeq2 library\n#library("DESeq2")')


# In[9]:


# Files to load into DE analysis (R)
metadata_file = str(paths.SAMPLE_ANNOT)
expression_data_file = str(paths.PAO1_GE_DE)
out_file = str(paths.DE_STATS)


# In[10]:


# Check ordering of sample ids
utils.check_sample_ordering(expression_data_file, metadata_file)


# In[11]:


"""%%R -i metadata_file -i expression_data_file -i out_file
expression_data <- t(as.matrix(read.csv(expression_data_file, sep="\t", header=TRUE, row.names=1)))
metadata <- as.matrix(read.csv(metadata_file, sep="\t", header=TRUE, row.names=1))

print("Checking sample ordering...")
print(all.equal(colnames(expression_data), rownames(metadata)))

group <- interaction(metadata[,1])

mm <- model.matrix(~0 + group)

#print(head(expression_data))

ddset <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~genotype)
print(ddset)
print(head( assay(ddset) ))
print(colData(ddset))
deseq_object <- DESeq(ddset)
print(deseq_object)

deseq_results <- results(deseq_object)
deseq_results

deseq_results_df <-  as.data.frame(deseq_results)

write.table(deseq_results_df, file = out_file, row.names = T, sep = "\t", quote = F)

# this is of class DESeqResults -- we want a data.frame
#deseq_df <- deseq_results %>%
  # make into data.frame
#  as.data.frame() %>%
  # the gene names are rownames -- let's make this it's own column for easy 
  # display
#  tibble::rownames_to_column(var = "Gene")

#deseq_df %>%
  # let's sort by statistic -- the highest values should be what is up in the
  # MYCN amplified cell lines
#  dplyr::arrange(dplyr::desc(stat))

#readr::write_tsv(deseq_df, path = deseq_df_file)"""


# In[13]:


get_ipython().run_cell_magic('R', '', 'library("limma")')


# In[14]:


get_ipython().run_cell_magic('R', '-i metadata_file -i expression_data_file -i out_file', "source('../core_acc_modules/DE_analysis.R')\n\nget_DE_stats(metadata_file,\n             expression_data_file,\n             out_file)")


# In[15]:


# Read in DE stats file
DE_stats = pd.read_csv(paths.DE_STATS, sep='\t', header=0, index_col=0)
print(DE_stats.shape)
DE_stats.head()


# ### Add gene names for ease of interpretation

# In[16]:


# Read gene number to name mapping
gene_name_mapping = pd.read_table(
    paths.GENE_ID2NAME,
    header=0,
    sep=',',
    index_col=0)

gene_name_mapping = gene_name_mapping[["Locus Tag", "Name"]]

gene_name_mapping.set_index("Locus Tag", inplace=True)
print(gene_name_mapping.shape)
gene_name_mapping.head()


# In[17]:


# Format gene numbers to remove extraneous quotes
gene_number = gene_name_mapping.index
gene_name_mapping.index = gene_number.str.strip("\"")

gene_name_mapping.dropna(inplace=True)
print(gene_name_mapping.shape)
gene_name_mapping.head(10)


# In[18]:


# Remove duplicate mapping
# Not sure which mapping is correct in this case
# PA4527 maps to pilC and still frameshift type 4 fimbrial biogenesis protein PilC (putative pseudogene)
gene_name_mapping = gene_name_mapping[~gene_name_mapping.index.duplicated(keep=False)]


# In[19]:


# Add gene names
#gene_name_mapping_dict = gene_name_mapping.to_dict()
DE_stats['Gene Name'] = DE_stats.index.map(gene_name_mapping["Name"])
DE_stats.head()


# In[20]:


# Plot
#DE_stats["-log10(padj)"] = -np.log10(DE_stats["padj"])
#sns.scatterplot(data=DE_stats, x="log2FoldChange", y="-log10(padj)")

DE_stats["-log10(padj)"] = -np.log10(DE_stats["adj.P.Val"])
sns.scatterplot(data=DE_stats, x="logFC", y="-log10(padj)")


# ### Compare out results with publication

# In[22]:


# Get number of DEGs
#selected_DE_stats = DE_stats[(DE_stats['padj']<0.05) & (abs(DE_stats['log2FoldChange'])>1)]
selected_DE_stats = DE_stats[(DE_stats['adj.P.Val']<0.05) & (abs(DE_stats['logFC'])>1)]

print(selected_DE_stats.shape)
selected_DE_stats.head()


# In[23]:


# Compare our findings against Sana et. al.
degs_sana = ["PA3431", "PA3432", "PA1244", "PA4685"]
for gene in degs_sana:
    if gene in list(selected_DE_stats.index):
        print(gene)


# In[24]:


# Compare our findings against Kim et. al.
'PA4236' in list(selected_DE_stats.index)


# **Conclusions:**
# * Our DE analysis found ~1K significantly differentially expressed genes
# * [Sana et. al](https://jb.asm.org/content/201/21/e00362-19) found ~ 2K DEGs between 2 strains where QS genes were DEGs. 
# * We found 2/4 genes highlighted in the Sana et. al. publication. But we're missing the main qsIA gene (PA1244) that the paper found to be more highly expressed in PAO1 vs PA14. One caveat to this publication is that its only looking at 2 strains so perhaps this explains the difference.
# * [Kim et. al.](https://link.springer.com/content/pdf/10.1007/s12275-019-9225-1.pdf) found katA (PA4236) is highly expressed in PAO1 but not PA14 strain. 
# * We did not find the katA gene in our analysis. One caveat to this publication is that its only looking at 2 strains so perhaps this explains the difference.
