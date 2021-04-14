"""Path definitions"""
from pathlib import Path

# Path to this repository
PROJECT_DIR = Path(__file__).parents[1]
ANALYSIS_DIR = PROJECT_DIR / "explore_data"
META_DIR = PROJECT_DIR / "data" / "metadata"

# Path to local directory where data files will be stored
LOCAL_DIR = Path.home()
LOCAL_DATA_DIR = LOCAL_DIR / "Documents" / "Data" / "Core_accessory" / "compendia_download"

# Location of gene expression matrix stored locally due to size
PAO1_GE = LOCAL_DATA_DIR / "TPM_PAO1_cdna_k15_cleaned.csv"
PA14_GE = LOCAL_DATA_DIR / "TPM_PA14_cdna_k15_cleaned.csv"

# Location of metadata mapping sample to strain name
SAMPLE_TO_STRAIN = LOCAL_DATA_DIR / "Run_Table_Strain_Bool_GD.csv"
SAMPLE_TO_STRAIN_PROCESSED = PROJECT_DIR / "data" / "metadata" / "SRA_annotations.tsv"

# Location for gene annotations from bactome, used to define core vs accessory genes
GENE_PAO1_ANNOT = LOCAL_DATA_DIR / "PAO1_ID_2_PA14_ID_PAO1ref.csv"
GENE_PA14_ANNOT = LOCAL_DATA_DIR / "PA14_ID_2_PAO1_ID_PA14ref.csv"

# Location for Salmon log files
PAO1_LOGS = LOCAL_DATA_DIR / "logs_pao1_cdna_k15.csv"
PA14_LOGS = LOCAL_DATA_DIR / "logs_pa14_cdna_k15.csv"

# Location for pre-binned compendia files
PAO1_PREBIN_COMPENDIUM = LOCAL_DATA_DIR / "pao1_prebin_compendia.tsv"
PA14_PREBIN_COMPENDIUM = LOCAL_DATA_DIR / "pa14_prebin_compendia.tsv"

# Location for processed compendia files
PAO1_COMPENDIUM_LABEL = LOCAL_DATA_DIR / "pao1_compendia_labeled.tsv"
PA14_COMPENDIUM_LABEL = LOCAL_DATA_DIR / "pa14_compendia_labeled.tsv"
PAO1_COMPENDIUM = LOCAL_DATA_DIR / "pao1_compendia.tsv"
PA14_COMPENDIUM = LOCAL_DATA_DIR / "pa14_compendia.tsv"

# Location of metadata
PAO1_REGULON = META_DIR / "regulons_format.csv"
PAO1_OPERON = META_DIR / "operons_format.csv"
