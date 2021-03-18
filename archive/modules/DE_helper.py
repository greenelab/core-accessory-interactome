"""
Author: Alexandra Lee
Date Created: 18 December 2020

Source code from
https://github.com/greenelab/generic-expression-patterns/blob/master/generic_expression_patterns_modules/DE_analysis.R

This script includes functions to prepare the data to run
DE analysis
"""

import pandas as pd


def process_samples_for_DESeq(
    expression_filename,
    grp_metadata_filename,
    out_expression_filename=None,
    count_threshold=None,
    process_metadata_filename=None,
):
    """
    This function processes samples in the template and simulated
    experiments to prepare for DE analysis using DESeq.

    These processing steps includes:
    1. Removing samples that are not included in the comparison.
    These "extra" samples occur when an experiment contains multiple
    comparisons.
    2. Removes genes with 0 counts across all samples
    3. (Optionally) filters genes with mean gene expression below
    some user defined threshold
    4. Case count values as integers
    5. Checks that the ordering of samples in the metadata file
    are consistent with the ordering in the gene expression data
    matrix. If the ordering is not consistent, then samples in
    the gene expression data matrix are re-ordered.

    Arguments
    ----------
    expression_filename: str
        File containing unnormalized gene expression data for
        either template or simulated experiments
    grp_metadata_filename: str
        File containing group assigments for samples to use
        for DESeq analysis
    out_expression_filename (optional): str
        File to save processed gene expression data to.
        If None then processed gene expression data will
        be output to the same input filename
    count_threshold (optinal): int
        Remove genes that have mean count <= count_threshold
        If None then no genes will be removed.
    process_metadata_filename (optional): str
        File containing assignment for which samples to drop.
        If None then all samples will be used.

    """

    # Read data
    expression = pd.read_csv(expression_filename, sep="\t", index_col=0, header=0)
    if process_metadata_filename is not None:
        process_metadata = pd.read_csv(
            process_metadata_filename, sep="\t", index_col=0, header=0
        )
    grp_metadata = pd.read_csv(grp_metadata_filename, sep="\t", header=0, index_col=0)

    if process_metadata_filename is not None:
        # Get samples ids to remove
        samples_to_remove = list(
            process_metadata[process_metadata["processing"] == "drop"].index
        )

        # Remove samples
        expression = expression.drop(samples_to_remove)

    # Cast as int
    expression = expression.astype(int)

    # Remove genes with 0 counts
    # all_zero_genes = list(expression.columns[(expression == 0).all()])
    # expression = expression.drop(columns=all_zero_genes)

    # assert len(list(expression.columns[(expression == 0).all()])) == 0

    # Remove genes below a certain threshold (if provided)
    if count_threshold is not None:
        genes_to_keep = expression.loc[:, expression.mean() >= count_threshold].columns
        expression = expression[genes_to_keep]

    # Check ordering of sample ids is consistent between gene expression data and metadata
    metadata_sample_ids = grp_metadata.index
    expression_sample_ids = expression.index

    if metadata_sample_ids.equals(expression_sample_ids):
        print("sample ids are ordered correctly")
    else:
        # Convert gene expression ordering to be the same as
        # metadata sample ordering
        print("sample ids don't match, going to re-order gene expression samples")
        expression = expression.reindex(metadata_sample_ids)

        assert expression.index.equals(metadata_sample_ids)

    # Save
    if out_expression_filename is not None:
        expression.to_csv(out_expression_filename, sep="\t")
    else:
        expression.to_csv(expression_filename, sep="\t")
