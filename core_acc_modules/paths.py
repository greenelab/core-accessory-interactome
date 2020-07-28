"""Path definitions"""
from pathlib import Path

# Path to this repository
PROJECT_DIR = Path(__file__).parents[1]
ANALYSIS_DIR = PROJECT_DIR / "sra_experiment"
METADATA_DIR = ANALYSIS_DIR / "data" / "metadata"
SAMPLE_ANNOT = METADATA_DIR / "sample_groups.txt"
GENE_PAO1_ANNOT = METADATA_DIR / "PAO1_ID_2_PA14_ID_PAO1ref.csv"
GENE_PA14_ANNOT = METADATA_DIR / "PA14_ID_2_PAO1_ID_PA14ref.csv"

# Path to local directory where data files will be stored
LOCAL_DIR = Path.home()

# Location where RNA-seq data is stored
NCBI_DIR = LOCAL_DIR / "ncbi" / "public"
SRA_DIR = NCBI_DIR / "sra"
FASTQ_DIR = NCBI_DIR / "fastq"

# Location where transcriptome references downloaded from Pseudomonas.com are stored
REF_DIR = LOCAL_DIR / "Documents" / "Data" / "Core_accessory"
PAO1_REF = REF_DIR / "Pseudomonas_aeruginosa_PAO1_107.ffn.gz"
PA14_REF = REF_DIR / "Pseudomonas_aeruginosa_UCBPP-PA14_109.ffn.gz"

# Location where mapping indices generated from `salmon index` are stored
PAO1_INDEX = REF_DIR / "pao1_index"
PA14_INDEX = REF_DIR / "pa14_index"

# Location where quantification results are stored from `salmon quant`
PAO1_QUANT = NCBI_DIR / "quants_pao1"
PA14_QUANT = NCBI_DIR / "quants_pa14"

# Location of gene expression matrix to use for correlation analysis
PAO1_GE = REF_DIR / "gene_expression_pao1_ref.tsv"
PA14_GE = REF_DIR / "gene_expression_pa14_ref.tsv"

# Location of gene expression matrix to use for DE analysis
PAO1_GE_DE = REF_DIR / "gene_expression_DE_input.tsv"
DE_STATS = REF_DIR / "DE_stats.tsv"

# Location of core genes to be reviewed
SHARED_CORE_PAO1_REF = REF_DIR / "shared_core_genes_pao1_ref.tsv"
PAO1_CORE_PAO1_REF = REF_DIR / "pao1_core_genes_pao1_ref.tsv"
PA14_CORE_PAO1_REF = REF_DIR / "pa14_core_genes_pao1_ref.tsv"

SHARED_CORE_PA14_REF = REF_DIR / "shared_core_genes_pa14_ref.tsv"
PAO1_CORE_PA14_REF = REF_DIR / "pao1_core_genes_pa14_ref.tsv"
PA14_CORE_PA14_REF = REF_DIR / "pa14_core_genes_pa14_ref.tsv"

# Location of accessory genes to be reviewed
PAO1_SAMPLE_PA14_REF = REF_DIR / "pao1_samples_pa14_ref.tsv"
PA14_SAMPLE_PAO1_REF = REF_DIR / "pa14_samples_pao1_ref.tsv"
HIGH_PA14_SAMPLE_PA14_REF = REF_DIR / "high_corr_pa14_samples_pa14_ref.tsv"
