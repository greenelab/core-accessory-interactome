"""Path definitions"""
from pathlib import Path

# Path to this repository
PROJECT_DIR = Path(__file__).parents[1]
ANALYSIS_DIR = PROJECT_DIR / "explore_data"
METADATA_DIR = ANALYSIS_DIR / "data" / "metadata"

# Path to local directory where data files will be stored
LOCAL_DIR = Path.home()
LOCAL_DATA_DIR = LOCAL_DIR / "Documents" / "Data" / "Core_accessory"

# Location of gene expression matrix to use for correlation analysis
PAO1_GE = LOCAL_DATA_DIR / "gene_expression_pao1_test_ref.tsv"
PA14_GE = LOCAL_DATA_DIR / "gene_expression_pa14_test_ref.tsv"

# Location of metadata
PAO1_METADATA = METADATA_DIR / "metadata_pao1.csv"
PA14_METADATA = METADATA_DIR / "metadata_pa14.csv"
