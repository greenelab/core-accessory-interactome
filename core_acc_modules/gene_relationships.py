"""
Author: Alexandra Lee
Date Created: 10 August 2021

Scripts to analyze relationships between gene types
"""
import pandas as pd
import numpy as np


def get_relationship_in_genome_space(core_acc_df, offset_to_bin, operon_df=None):
    """
    For each accessory gene, is the 1-NN/2-NN/3-NN core or accessory? This function
    calculates the number of accessory-accessory relationships that are 1-away, 2-away, etc.
    This function also calculates the number of accessory-core relationships that are
    1-away, 2-away, etc. These same calculations are performed starting with a core gene.

    Arguments:
    ----------
    core_acc_df: pandas df
        Dataframe with gene ids on the index and core/acc label as a column
    offset_to_bin: int
        Calculations will increment by 1 up until offset_to_bin. Any counts corresponding
        to nearest neighbors that exceed the `offset_to_bin` will be averaged.
    operon_df: pandas df
        Dataframe with gene ids on the index and operon name as a column.
        If this argument is provided, then genes within the same operon
        will not be included in the counts.
    """
    gene_type_start = ["acc", "core"]
    gene_type_compare = ["acc", "core"]

    core_acc_df_len = len(core_acc_df)
    offset_max = core_acc_df_len - 1

    core_acc_df_pad = np.pad(
        core_acc_df["core/acc"], offset_max, "constant", constant_values="NA"
        )
    if operon_df is not None:
        operon_df_pad = np.pad(
            operon_df["operon_name"], offset_max, "constant", constant_values="NA"
        )

    assert core_acc_df.shape == operon_df.shape

    rows = []
    for gene_start in gene_type_start:
        for gene_compare in gene_type_compare:
            for offset in range(1, offset_max + 1):
                # Print statements to understand what is happening
                # print("start left", core_acc_df_pad[offset_max : core_acc_df_len + offset_max])
                # print("compare left", core_acc_df_pad[
                #             offset_max - offset : core_acc_df_len + offset_max - offset
                #         ])
                #
                # print("start left", operon_df_pad[offset_max : core_acc_df_len + offset_max])
                # print("compare left", operon_df_pad[
                #             offset_max - offset : core_acc_df_len + offset_max - offset
                #         ])

                # Compare left nearest neighbors: Are they core or accessory?
                core_acc_left = (
                    core_acc_df_pad[offset_max : core_acc_df_len + offset_max]
                    == gene_start
                    ) & (
                        core_acc_df_pad[
                            offset_max - offset : core_acc_df_len + offset_max - offset
                            ]
                            == gene_compare
                            )

                # Compare right nearest neighbors: Are they core or accessory?
                core_acc_right = (
                    core_acc_df_pad[offset_max: core_acc_df_len + offset_max]
                    == gene_start
                    ) & (
                        core_acc_df_pad[
                            offset_max + offset: core_acc_df_len + offset_max + offset
                            ]
                            == gene_compare
                            )

                if operon_df is not None:
                    # Compare left operons: Are the genes in the same operon?
                    operon_left = (
                        operon_df_pad[offset_max: core_acc_df_len + offset_max]
                        == operon_df_pad[
                            offset_max - offset: core_acc_df_len + offset_max - offset
                        ]
                    )

                    # Compare right operons: Are the genes in the same operon?
                    operon_right = (
                        operon_df_pad[offset_max: core_acc_df_len + offset_max]
                        == operon_df_pad[
                            offset_max + offset: core_acc_df_len + offset_max + offset
                        ]
                    )

                    # Sum all comparisons
                    counts = (core_acc_left & ~operon_left).sum() + (
                        core_acc_right & ~operon_right
                    ).sum()
                else:
                    counts = core_acc_left.sum() + core_acc_right.sum()

                rows.append(
                    {
                        "gene start": gene_start,
                        "gene compare": gene_compare,
                        "offset": offset,
                        "total": counts,
                    }
                )

    genome_dist_counts = pd.DataFrame(rows)

    # Bin distances above offset_to_bin
    long_dist = (
        genome_dist_counts.query("offset>@offset_to_bin")
        .groupby(["gene compare", "gene start"])["total"]
        .mean()
        .reset_index()
    )
    long_dist["offset"] = f"{offset_to_bin}+"
    genome_dist_counts = genome_dist_counts.query("offset<=@offset_to_bin").append(
        long_dist, ignore_index=True
    )

    return genome_dist_counts


def get_relationship_in_expression_space(
    corr_df,
    genes_to_consider,
    gene_mapping_df,
    offset_to_bin,
    operon_df=None,
    sum_increment=1,
):

    """For each accessory gene, is the highest correlated, 2nd-highest correlated, 3rd
    highest correlated gene, etc core or accessory? This function
    calculates the number of accessory-accessory relationships that are highest correlated,
    next highest correlated, third highest correlated, etc.
    This function also calculates the number of accessory-core relationships that are
    highest correlated, next highest correlated, third highest correlated, etc.
    These same calculations are performed starting with a core gene.

    Arguments:
    ----------
    corr_df: pandas df
        Dataframe containing the correlation matrix. The matrix is a gene x gene
        dataframe containing values that correspond to how correlated each pair of
        genes is
    genes_to_consider: list
        List of gene ids that will be used as starting genes
    gene_mapping_df: pandas df
        Dataframe with gene ids on the index and core/acc label as a column
    offset_to_bin: int
        Calculations will increment by 1 up until offset_to_bin. Any counts corresponding
        to nearest neighbors that exceed the `offset_to_bin` will be averaged.
    operon_df: pandas df
        Dataframe with gene ids on the index and operon name as a column.
        If this argument is provided, then genes within the same operon
        will not be included in the counts.
    sum_increment: int
        Integer that is the window size to sum across. This will allow
        users to consider the top N most correlated genes instead of just
        the top correlated gene.
    """
    # Get subset of genes
    corr_subset = corr_df.loc[genes_to_consider]

    # Note: PA14 contains duplicate rows so we will drop those here
    corr_subset = corr_subset.drop_duplicates()

    rows = []
    for gene in corr_subset.index:
        if operon_df is not None:
            # This subset needs to be reset each iteration
            # since we are dropping columns below
            corr_subset = corr_df.loc[genes_to_consider]

            # Note: PA14 contains duplicate rows so we will drop those here
            corr_subset = corr_subset.drop_duplicates()

            # Check if gene is found in an operon
            if gene in operon_df.index:
                # Find operons containing 'gene'
                group_name = operon_df.loc[gene, "operon_name"]

                # Dictionary format: pao1_operon_dict[operon_name] = [list of genes]
                operon_dict = operon_df.groupby("operon_name").groups
                co_operonic_genes = list(operon_dict[group_name])

                # Remove columns corresponding to co-operonic genes from the corr_subset
                co_operonic_genes_to_remove = list(
                    set(corr_subset.columns).intersection(co_operonic_genes)
                )
                corr_subset = corr_subset.drop(columns=co_operonic_genes_to_remove)

        offset_max = corr_subset.shape[1]

        top_corr_genes = list(corr_subset.loc[gene].nlargest(offset_max).index[1:])
        top_gene_labels = list(gene_mapping_df.loc[top_corr_genes, "core/acc"].values)
        rows.append(top_gene_labels)

    expression_dist_counts = pd.DataFrame(rows)

    # Count types of relationships
    expression_dist_counts_acc = (expression_dist_counts == "acc").sum().to_frame("acc")

    expression_dist_counts = expression_dist_counts_acc.join(
        (expression_dist_counts == "core").sum().to_frame("core")
    )

    if sum_increment > 1:
        expression_dist_counts = (
            expression_dist_counts.rolling(sum_increment)
            .sum()
            .iloc[::-1]
            .shift(1)
            .sort_index()
        )

    # Format counts for plotting
    expression_dist_counts = expression_dist_counts.melt(
        var_name="gene type", value_name="total", ignore_index=False
    )
    expression_dist_counts = expression_dist_counts.rename_axis("offset").reset_index()
    expression_dist_counts["offset"] = expression_dist_counts["offset"] + 1

    # Average counts for weaker correlation relationships
    weak_corr = (
        expression_dist_counts.query("offset>@offset_to_bin")
        .groupby("gene type")["total"]
        .mean()
        .to_frame()
    )
    weak_corr = weak_corr.reset_index()
    weak_corr["offset"] = f"+{offset_to_bin}"

    expression_dist_counts = expression_dist_counts.query(
        "offset<=@offset_to_bin"
        ).append(weak_corr, ignore_index=True)

    # Add proportion - How should we calculate proportion
    # Of all 1-NN, %accessory, %core
    # Or, of all the accessory genes, %are 1-NN, of all core genes, %are 1-NN
    # total_counts = expression_dist_counts.groupby("offset")["total"].sum()[1]
    # expression_dist_counts["proportion"] = expression_dist_counts["total"]/total_counts
    return expression_dist_counts
