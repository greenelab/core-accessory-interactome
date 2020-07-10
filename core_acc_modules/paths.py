"""Path definitions"""
from pathlib import Path

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
