"""Path definitions"""
from pathlib import Path

# Path to this repository
PROJECT_DIR = Path(__file__).parents[1]
ANALYSIS_DIR = PROJECT_DIR / "explore_data"

# Path to local directory where data files will be stored
LOCAL_DIR = Path.home()
LOCAL_DATA_DIR = LOCAL_DIR / "Documents" / "Data" / "Core_accessory" / "compendia_download"

# Location of gene expression matrix stored locally due to size
PAO1_GE = LOCAL_DATA_DIR / "rnaseq_compendium_cleaned.csv"
PA14_GE = LOCAL_DATA_DIR / "TPM_pa14_cdna_k15.csv"

# Location of metadata mapping sample to strain name
SAMPLE_TO_STRAIN = LOCAL_DATA_DIR / "Run_Table_Strain_Bool_GD.csv"

# Location for gene annotations from bactome, used to define core vs accessory genes
GENE_PAO1_ANNOT = LOCAL_DATA_DIR / "PAO1_ID_2_PA14_ID_PAO1ref.csv"
GENE_PA14_ANNOT = LOCAL_DATA_DIR / "PA14_ID_2_PAO1_ID_PA14ref.csv"
