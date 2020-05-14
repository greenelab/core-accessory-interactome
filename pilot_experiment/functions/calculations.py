'''
Author: Alexandra Lee
Date Created: 17 April 2020

Scripts to process expression data for pilot experiment
'''
import pandas as pd
import os
import random
import numpy as np
import pickle

import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
    
np.random.seed(123)

def get_coexpression_stats(correlation_file,
operon_annotation_file,
core_gene_ids_file,
acc_gene_ids_file,
coexpression_threshold):
    """
    Returns summary statistics about co-expressed and co-operonic genes

    Arguments
    ----------
    correlation_file: str
        Pickle file containing pairwise correlation between all genes

    operon_annotation_file: str
        File containing list of genes per operon

    core_gene_id_file: str
        File containing core genes annotated from annotate_genes() in process_data.py

    acc_gene_id_file: str
        File containing accessory genes annotated from annotate_genes() in process_data.py

    coexpression_threshold: float
        Hard threshold cutoff that determine if genes are co-expressed .
        If corr(gene i, gene j) > coexpression_threshold, then gene i and j are co-expressed 
    """
    # Read correlation data
    core_gene_ids = pickle.load(open(core_gene_ids_file, "rb"))
    acc_gene_ids = pickle.load(open(acc_gene_ids_file, "rb"))
    corr_data = pickle.load(open(correlation_file, "rb"))

    # Get all gene ids
    all_gene_ids = list(corr_data.index)

    # Define threshold for highly co-expressed genes
    threshold = coexpression_threshold

    # Apply threshold to identify which genes are co-expressed
    coexpressed = corr_data>threshold

    # Get total number of genes that are co-expressed per gene
    # This will serve as a test case
    num_coexpressed_genes = coexpressed.sum(axis=1)

    # Given the list of co-expressed genes
    # we want to differentiate between those that are core and those that are accessory
    name_coexpressed_genes = {}
    num_coexpressed_core_genes = {}
    num_coexpressed_acc_genes = {}

    for gene_id in all_gene_ids:
        # Get row of correlation matrix
        # The values in the row corresponds to if there exists a gene is co-expressed with the gene_id
        coexpressed_gene_values = coexpressed.loc[gene_id]
        
        # Check that our calculations are consistent
        assert(num_coexpressed_genes[gene_id] == sum(coexpressed_gene_values))
        
        if num_coexpressed_genes[gene_id] > 0:
            # Get list of co-expressed genes
            lst_coexpressed_genes = list(coexpressed_gene_values[coexpressed_gene_values].index)
            name_coexpressed_genes[gene_id] = lst_coexpressed_genes
            
            # Get the number of co-expressed genes that are core, accessory
            num_core_genes = len(set(lst_coexpressed_genes).intersection(core_gene_ids))
            num_acc_genes = len(set(lst_coexpressed_genes).intersection(acc_gene_ids))
            num_coexpressed_core_genes[gene_id] = num_core_genes
            num_coexpressed_acc_genes[gene_id] = num_acc_genes
            
        else:
            name_coexpressed_genes[gene_id] = []
            num_coexpressed_core_genes[gene_id] = 0
            num_coexpressed_acc_genes[gene_id] = 0
    
    # Calculate ratio of core:accessory genes in co-expressed gene sets
    coexpressed_core_prop = {}
    coexpressed_acc_prop = {}

    for gene_id in all_gene_ids:
        num_core_genes = num_coexpressed_core_genes[gene_id]
        num_acc_genes = num_coexpressed_acc_genes[gene_id]
        if (num_core_genes == 0 & num_acc_genes == 0):
            coexpressed_core_prop[gene_id] = 0
            coexpressed_acc_prop[gene_id] = 0
        else:
            coexpressed_core_prop[gene_id] = num_core_genes/(num_core_genes + num_acc_genes)
            coexpressed_acc_prop[gene_id] = num_acc_genes/(num_core_genes + num_acc_genes)

    # Get operon annotations per gene id
    name_cooperonic_genes = map_gene_to_cooperonic_genes(operon_annotation_file,
    all_gene_ids)

    # Compare co-expressed gene set and co-operonic genes per reference gene id
    num_non_cooperonic_coexpressed_genes = {}
    num_non_cooperonic_coexpressed_core_genes = {}
    num_non_cooperonic_coexpressed_acc_genes = {}

    for gene_id in all_gene_ids:
        # Get co-operonic gene list
        cooperonic_genes = name_cooperonic_genes[gene_id]
        
        # Get co-expressed gene list
        coexpressed_genes = name_coexpressed_genes[gene_id]
        
        # Find non co-operonic genes
        # Find genes that DO NOT intersect between co-operonic genes and co-expressed genes
        cooperonic_coexpressed_genes = set(coexpressed_genes).intersection(cooperonic_genes)
        
        non_cooperonic_coexpressed_genes = set(coexpressed_genes) - cooperonic_coexpressed_genes
        
        # Get number of non-co-operonic genes
        num_non_cooperonic_coexpressed_genes[gene_id] = len(non_cooperonic_coexpressed_genes)
        
        if num_non_cooperonic_coexpressed_genes[gene_id] > 0:        
            # Get the number of non co-operonic co-expressed genes that are core, accessory
            num_core_genes = len(non_cooperonic_coexpressed_genes.intersection(core_gene_ids))
            num_acc_genes = len(non_cooperonic_coexpressed_genes.intersection(acc_gene_ids))
            num_non_cooperonic_coexpressed_core_genes[gene_id] = num_core_genes
            num_non_cooperonic_coexpressed_acc_genes[gene_id] = num_acc_genes
            
        else:
            num_non_cooperonic_coexpressed_core_genes[gene_id] = 0
            num_non_cooperonic_coexpressed_acc_genes[gene_id] = 0
    
    # Calculate ratio of core:accessory genes in co-expressed gene sets
    non_cooperonic_coexpressed_core_prop = {}
    non_cooperonic_coexpressed_acc_prop = {}

    for gene_id in all_gene_ids:
        num_core_genes = num_non_cooperonic_coexpressed_core_genes[gene_id]
        num_acc_genes = num_non_cooperonic_coexpressed_acc_genes[gene_id]
        if (num_core_genes == 0 & num_acc_genes == 0):
            non_cooperonic_coexpressed_core_prop[gene_id] = 0
            non_cooperonic_coexpressed_acc_prop[gene_id] = 0
        else:
            non_cooperonic_coexpressed_core_prop[gene_id] = num_core_genes/(num_core_genes + num_acc_genes)
            non_cooperonic_coexpressed_acc_prop[gene_id] = num_acc_genes/(num_core_genes + num_acc_genes)

    # Core gene stats
    core_stats_df = pd.DataFrame(data={'ref_gene':core_gene_ids,
                                    'num_coexpressed_genes':num_coexpressed_genes[core_gene_ids],
                                    'num_coexpressed_core': [num_coexpressed_core_genes[k] for k in core_gene_ids],
                                    'num_coexpressed_acc': [num_coexpressed_acc_genes[k] for k in core_gene_ids],
                                    'percent_coexpressed_core': [coexpressed_core_prop[k] for k in core_gene_ids],
                                    'percent_coexpressed_acc': [coexpressed_acc_prop[k] for k in core_gene_ids],
                                    'num_non_cooperonic_coexpressed_genes':[num_non_cooperonic_coexpressed_genes[k] 
                                                                            for k in core_gene_ids],
                                    'num_non_cooperonic_coexpressed_core': [num_non_cooperonic_coexpressed_core_genes[k] 
                                                                            for k in core_gene_ids],
                                    'num_non_cooperonic_coexpressed_acc': [num_non_cooperonic_coexpressed_acc_genes[k] 
                                                                            for k in core_gene_ids],
                                    'percent_non_cooperonic_coexpressed_core': [non_cooperonic_coexpressed_core_prop[k] 
                                                                                for k in core_gene_ids],
                                    'percent_non_cooperonic_coexpressed_acc': [non_cooperonic_coexpressed_acc_prop[k] 
                                                                                for k in core_gene_ids]
                                    }
                                )

    # Accessory gene stats
    acc_stats_df = pd.DataFrame(data={'ref_gene':acc_gene_ids,
                                    'num_coexpressed_genes':num_coexpressed_genes[acc_gene_ids],
                                    'num_coexpressed_core': [num_coexpressed_core_genes[a] for a in acc_gene_ids],
                                    'num_coexpressed_acc': [num_coexpressed_acc_genes[a] for a in acc_gene_ids],
                                    'percent_coexpressed_core': [coexpressed_core_prop[a] for a in acc_gene_ids],
                                    'percent_coexpressed_acc': [coexpressed_acc_prop[a] for a in acc_gene_ids],
                                    'num_non_cooperonic_coexpressed_genes':[num_non_cooperonic_coexpressed_genes[a] 
                                                                            for a in acc_gene_ids],
                                    'num_non_cooperonic_coexpressed_core': [num_non_cooperonic_coexpressed_core_genes[a] 
                                                                            for a in acc_gene_ids],
                                    'num_non_cooperonic_coexpressed_acc': [num_non_cooperonic_coexpressed_acc_genes[a] 
                                                                            for a in acc_gene_ids],
                                    'percent_non_cooperonic_coexpressed_core': [non_cooperonic_coexpressed_core_prop[a] 
                                                                                for a in acc_gene_ids],
                                    'percent_non_cooperonic_coexpressed_acc': [non_cooperonic_coexpressed_acc_prop[a] 
                                                                                for a in acc_gene_ids]
                                    }
                                )
    return core_stats_df, acc_stats_df
            

def map_gene_to_cooperonic_genes(operon_annotation_file,
all_gene_ids):
    """
    Returns mapping between gene id and list of co-operonic genes

    Arguments
    ----------
    operon_annotation_file: str
        File containing list of genes per operon

    all_gene_ids: list
        List of all gene ids to search for co-operonic genes 
    """

    # Read operon data
    # Manually had to set names to be the max size operon
    operon_data = pd.read_csv(
        operon_annotation_file,
        header=None,
        sep='\t',
        names=range(15)
    )
    # Associate operons with reference gene
    name_cooperonic_genes = {}

    for gene_id in all_gene_ids:
        # Get operons containing reference gene_id   
        # Search for gene_id in each operon
        operon_search = operon_data.where(operon_data == gene_id).dropna(how='all').dropna(axis=1)
        
        # Note: this annotation file has a 1:1 mappring between gene and list of co-operonic genes
        if operon_search.empty:
            name_cooperonic_genes[gene_id] = []
        else:
            row_id = operon_search.index[0]
            name_cooperonic_genes[gene_id] = list(operon_data.loc[row_id].dropna())

    return name_cooperonic_genes