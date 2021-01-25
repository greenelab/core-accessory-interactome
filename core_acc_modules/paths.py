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

# Location where raw transcriptome references downloaded from Pseudomonas.com andn NCBI are stored
REF_DIR = LOCAL_DIR / "Documents" / "Data" / "Core_accessory"
PAO1_REF = REF_DIR / "Pseudomonas_aeruginosa_PAO1_107.fasta"
PA14_REF = REF_DIR / "Pseudomonas_aeruginosa_UCBPP-PA14_109.fasta"
PHAGE_REF = REF_DIR / "phage_sequences.fasta"

# Location of BLAST DB
BLAST_DIR = REF_DIR / "blast" / "db"
PAO1_DB_DIR = BLAST_DIR / "PAO1_DB"
PA14_DB_DIR = BLAST_DIR / "PA14_DB"
PAO1_BLAST_RESULT = BLAST_DIR / "pao1_blast_output.tsv"
PA14_BLAST_RESULT = BLAST_DIR / "pa14_blast_output.tsv"

# Location processed references
PAO1_PHAGE_REF = REF_DIR / "Pseudomonas_aeruginosa_PAO1_107_phage.fasta"
PA14_PHAGE_REF = REF_DIR / "Pseudomonas_aeruginosa_UCBPP-PA14_109_phage.fasta"

# Location where mapping indices generated from `salmon index` are stored
PAO1_PHAGE_INDEX = REF_DIR / "pao1_phage_index"
PA14_PHAGE_INDEX = REF_DIR / "pa14_phage_index"

# Location where quantification results are stored from `salmon quant`
PAO1_PHAGE_QUANT = NCBI_DIR / "quants_pao1_phage"
PA14_PHAGE_QUANT = NCBI_DIR / "quants_pa14_phage"

# Location of gene expression matrix to use for correlation analysis
PAO1_PHAGE_GE = REF_DIR / "gene_expression_pao1_phage_ref.tsv"
PA14_PHAGE_GE = REF_DIR / "gene_expression_pa14_phage_ref.tsv"
