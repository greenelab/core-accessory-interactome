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
PAO1_GE = LOCAL_DATA_DIR / "pao1_aligned_rnaseq_compendium_zp2_MRnorm.csv"
PA14_GE = LOCAL_DATA_DIR / "pa14_aligned_rnaseq_compendium_zp2_MRnorm.csv"

# Location of metadata mapping sample to strain name
SAMPLE_TO_STRAIN = LOCAL_DATA_DIR / "Run_Table_Strain_Bool_GD.csv"
SAMPLE_TO_STRAIN_PROCESSED = PROJECT_DIR / "data" / "metadata" / "SRA_annotations.tsv"

# Location of metadata for samples including biosample, instrument, center
SAMPLE_METADATA = LOCAL_DATA_DIR / "SraRunTable.csv"

# Location for gene annotations from bactome, used to define core vs accessory genes
GENE_PAO1_ANNOT = LOCAL_DATA_DIR / "PAO1_ID_2_PA14_ID_PAO1ref.csv"
GENE_PA14_ANNOT = LOCAL_DATA_DIR / "PA14_ID_2_PAO1_ID_PA14ref.csv"

# Location for Salmon log files
PAO1_LOGS = LOCAL_DATA_DIR / "logs_pao1_cdna_k15.csv"
PA14_LOGS = LOCAL_DATA_DIR / "logs_pa14_cdna_k15.csv"

# Location for pre-binned compendia files that have been formated to be
# sample x gene matrices with experiment id as sample id
PAO1_PREBIN_COMPENDIUM = LOCAL_DATA_DIR / "pao1_prebin_compendia.tsv"
PA14_PREBIN_COMPENDIUM = LOCAL_DATA_DIR / "pa14_prebin_compendia.tsv"

# Location for processed/binned compendia files
PAO1_COMPENDIUM_LABEL = LOCAL_DATA_DIR / "pao1_compendia_labeled.tsv"
PA14_COMPENDIUM_LABEL = LOCAL_DATA_DIR / "pa14_compendia_labeled.tsv"
PAO1_COMPENDIUM = LOCAL_DATA_DIR / "pao1_compendia.tsv"
PA14_COMPENDIUM = LOCAL_DATA_DIR / "pa14_compendia.tsv"

# Location of correlation matrices
PAO1_CORR_RAW = LOCAL_DATA_DIR / "pao1_pearson_mat.tsv"
PA14_CORR_RAW = LOCAL_DATA_DIR / "pa14_pearson_mat.tsv"
PAO1_CORR_LOG_SPELL = LOCAL_DATA_DIR / "pao1_all_log_spell_mat.tsv"
PA14_CORR_LOG_SPELL = LOCAL_DATA_DIR / "pa14_all_log_spell_mat.tsv"
PAO1_CORR_LOG_SPELL_CORE = LOCAL_DATA_DIR / "pao1_core_log_spell_mat.tsv"
PA14_CORR_LOG_SPELL_CORE = LOCAL_DATA_DIR / "pa14_core_log_spell_mat.tsv"
PAO1_CORR_LOG_SPELL_ACC = LOCAL_DATA_DIR / "pao1_acc_log_spell_mat.tsv"
PA14_CORR_LOG_SPELL_ACC = LOCAL_DATA_DIR / "pa14_acc_log_spell_mat.tsv"

# Location of metadata
PAO1_REGULON = META_DIR / "regulons_format.csv"
PAO1_OPERON = META_DIR / "PAO1-operons-2021-07-19.csv"
PA14_OPERON = META_DIR / "PA14-operons-2021-07-19.csv"

# Location of PAO1 array compendium metadata file
ARRAY_DATA_URL = "https://raw.githubusercontent.com/greenelab/adage/master/Data_collection_processing/Pa_compendium_02.22.2014.pcl"
ARRAY_METADATA_URL = "https://raw.githubusercontent.com/greenelab/generic-expression-patterns/97a55c8d53b5d1812399479d530b9cbaee689079/pseudomonas_analysis/data/metadata/sample_annotations.tsv"

# Location of processed PAO1 array and RNA-seq compendia files
# These compendia are using the same set of genes to compare
# module composition between array and RNA-seq
ARRAY_COMPENDIUM_TO_COMPARE = LOCAL_DATA_DIR / "pao1_array_compendia_tocompare.tsv"
RNASEQ_COMPENDIUM_TO_COMPARE = LOCAL_DATA_DIR / "pao1_rnaseq_compendia_tocompare.tsv"

# Location of SOPHIE identified common DEGs
COMMON_DEGS_PAO1 = LOCAL_DATA_DIR / "generic_gene_summary_SRP117105.tsv"
COMMON_DEGS_PA14 = LOCAL_DATA_DIR / "generic_gene_summary_SRP074292.tsv"
