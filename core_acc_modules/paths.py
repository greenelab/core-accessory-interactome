"""Path definitions"""
from pathlib import Path

# Path to this repository
PROJECT_DIR = Path(__file__).parents[1]
ANALYSIS_DIR = PROJECT_DIR / "test_phage_reference"
METADATA_DIR = ANALYSIS_DIR / "data" / "metadata"
SRA_ACC = METADATA_DIR / "sra_acc.txt"

# Path to local directory where data files will be stored
LOCAL_DIR = Path.home()

# Location where RNA-seq data is stored
NCBI_DIR = LOCAL_DIR / "ncbi" / "public"
SRA_DIR = NCBI_DIR / "sra"
FASTQ_DIR = NCBI_DIR / "fastq_phage"

# Location where transcriptome references downloaded from Pseudomonas.com are stored
REF_DIR = LOCAL_DIR / "Documents" / "Data" / "Core_accessory"
RAW_PHAGE_REF = REF_DIR / "phage_sequences.fasta"
PAO1_REF = REF_DIR / "Pseudomonas_aeruginosa_PAO1_107.ffn.gz"
PA14_REF = REF_DIR / "Pseudomonas_aeruginosa_UCBPP-PA14_109.ffn.gz"
PHAGE_REF = REF_DIR / "phage_sequences_processed.fasta"

# Location where mapping indices generated from `salmon index` are stored
PAO1_INDEX = REF_DIR / "pao1_index"
PA14_INDEX = REF_DIR / "pa14_index"
PHAGE_INDEX = REF_DIR / "phage_index"

# Location where quantification results are stored from `salmon quant`
PAO1_QUANT = NCBI_DIR / "quants_pao1"
PA14_QUANT = NCBI_DIR / "quants_pa14"
PHAGE_QUANT = NCBI_DIR / "quants_phage"

# Location of gene expression matrix to use for correlation analysis
PAO1_GE = REF_DIR / "gene_expression_pao1_ref.tsv"
PA14_GE = REF_DIR / "gene_expression_pa14_ref.tsv"
PHAGE_GE = REF_DIR / "gene_expression_phage_ref.tsv"
