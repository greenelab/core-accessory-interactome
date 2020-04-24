'''
Author: Alexandra Lee
Date Created: 17 April 2020

Scripts to process expression data for pilot experiment
'''
import pandas as pd
import os
import random

import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

def select_expression_data(input_data_file,
metadata_file,
lst_experiments,
output_file):
    """
    Returns a subset of expression data, specified by the experiment_ids
    provided by the user, to use for downstream exploration analysis

    Arguments
    ----------
    input_data_file: str
        File containing normalized gene expression data

        ------------------------------| PA0001 | PA0002 |...
        05_PA14000-4-2_5-10-07_S2.CEL | 0.8533 | 0.7252 |...
        54375-4-05.CEL                | 0.7789 | 0.7678 |...
        ...                           | ...    | ...    |...
    
    metadata_file: str
        File containing the mapping between experiment ids and sample ids. This 
        mapping can be used to link sample ids in the input_data_file with the 
        experiment ids provided by the user

    lst_experiments: list
        List of experiment ids

    output_file: str
        Filename to save selected expression data

    """
    # Read data
    normalized_data = pd.read_csv(
        input_data_file,
        header=0,
        sep='\t',
        index_col=0).T

    # Read metadata
    metadata = pd.read_csv(
        metadata_file,
        header=0,
        sep='\t',
        index_col=0)
        
    map_experiment_sample = metadata[['sample_name', 'ml_data_source']]
    
    # Get expression data associated with experiment_id
    selected_mapping = map_experiment_sample.loc[lst_experiments]
    selected_sample_ids = list(selected_mapping['ml_data_source'].values)

    selected_data = normalized_data.loc[selected_sample_ids]

    print('The selected dataset contains {} samples and {} genes'.
          format(selected_data.shape[0],selected_data.shape[1]))

    # Save selected gene expression data
    selected_data.to_csv(output_file, sep='\t', index=True)

def permute_expression_data(input_data_file,
output_file):
    """
    Returns permuted version of input expression data to be used as a baseline
    for downstream analysis

    Arguments
    ----------
    input_data_file: str
        File containing normalized gene expression data

        ------------------------------| PA0001 | PA0002 |...
        05_PA14000-4-2_5-10-07_S2.CEL | 0.8533 | 0.7252 |...
        54375-4-05.CEL                | 0.7789 | 0.7678 |...
        ...                           | ...    | ...    |...

    output_file: str
        Filename to save permuted expression data
    """
    # Read data
    selected_data = pd.read_csv(
        input_data_file,
        header=0,
        sep='\t',
        index_col=0)

    # For each sample, shuffle gene values (i.e. shuffle values independently per sample) in order to
    # disrupt any gene structure
    shuffled_arr = []
    num_samples = selected_data.shape[0]

    for i in range(num_samples):
        row = list(selected_data.values[i])
        shuffled_row = random.sample(row, len(row))
        shuffled_arr.append(shuffled_row)

    shuffled_selected_data = pd.DataFrame(shuffled_arr,
                                        index=selected_data.index,
                                        columns=selected_data.columns)

    # Save permuted gene expression data
    shuffled_selected_data.to_csv(output_file, sep='\t', index=True)

def annotate_genes(input_data_file,
gene_annotation_file,
output_file):
    """
    Returns file mapping gene ids to label ('core' or 'accessory') using gene annotations from
    BACTOME annotation tool: https://pseudomonas-annotator.shinyapps.io/pa_annotator/

     Arguments
    ----------
    input_data_file: str
        File containing normalized gene expression data

        ------------------------------| PA0001 | PA0002 |...
        05_PA14000-4-2_5-10-07_S2.CEL | 0.8533 | 0.7252 |...
        54375-4-05.CEL                | 0.7789 | 0.7678 |...
        ...                           | ...    | ...    |...

    gene_annotation_file: str
        File containing mapping between PAO1 and PA14 gene ids

        Columns: PAO1_ID,Name,Product.Name,PA14_ID

    output_file: str
        Filename to save gene mapping

    """
    # Read data
    selected_data = pd.read_csv(
        input_data_file,
        header=0,
        sep='\t',
        index_col=0)

    all_genes = list(selected_data.columns)

    # Read gene annotation 
    gene_mapping = pd.read_csv(
        gene_annotation_file)

    # Accessory are genes that don't have a mapping to the PA14 ID
    unmapped_genes = gene_mapping[gene_mapping["PA14_ID"].isna()]

    PAO1_only = list(unmapped_genes["PAO1_ID"])
    print("No. of PAO1 only genes: {}".format(len(PAO1_only)))

    # Create df with PAO1 gene IDs and label for core/accessory
    annot =  ["accessory" if (gene_id in PAO1_only) else "core" for gene_id in all_genes]

    df = pd.DataFrame(list(zip(all_genes, annot)), 
                columns =['PAO1_gene_id', 'annotation']) 

    # Save labels
    df.to_csv(output_file, sep='\t', index=False)
