# Compare SPELL vs counts correlation

There are two possible correlation matrices we can use for our analysis:

1. Correlation of the MR counts expression matrix
2. Correlation of the MR counts that have been processed using SPELL

**About the correlation matrices**

The correlation of the counts matrix relates how similar a pair of genes are based on their expression profiles - **relates genes over samples**.
High correlation means that a pair of genes have a similar expression profiles - i.e. similar levels of expression across samples/contexts, so genes both have low expression in the same samples and high expression in the same samples.
* Pro: Easy to interpret
* Con: Many gene pairs found to have a high correlation because many genes are related to the same pathway have the similar expression profiles. This is consistent with [Myers et al.](https://link.springer.com/article/10.1186/1471-2164-7-187), who found that there can be an over-representation of genes associated with the same pathway (i.e. a large fraction of gene pairs represent ribosomal relationships). This very prominent signal makes it difficult to detect other signals. Figure 1C demonstrates that a large fraction of gene pairs are ribosomal relationships - in the top 0.1% most co-expressed genes, 99% belong to the ribosome pathway. Furthermore, protein function prediction based on co-expression drop dramatically after removing the ribisome pathway (Figure 1A, B).

To try to remove this very dominant global signal in the data. Here we are applying dimensionality reduction techniques in addition to scaling the data using a method called SPELL.
The correlation of the SPELL matrix relates genes based on the gene coefficient matrix - **relate genes over their contribution to singular vectors (linear combination of genes - linear relationship between genes)**.
High correlation means that a pair of genes contributes similarly to a singular vector, which are the axes pointing in the direction of the spread of the data and capture how genes are related to each other
* Pro: Gene contributions are more balanced so that redundant signals (i.e. many genes from the same pathway - genes that vary together) are represented by a few SVs as opposed to many samples. More balanced also means that more subtle signals can be amplified (i.e. genes related by a smaller pathway are also captured by a few SVs)
* Con: Can amplify noise - i.e. an SV that corresponds to some technical source of variability now has a similar weight to other real signals

For more information comparing using counts vs SPELL-processing see: https://docs.google.com/presentation/d/18E0boNODJaxP-YYNIlccrh0kASbc7bapQBMovOX62jw/edit#slide=id.gf9d09c6be6_0_0

These notebooks generate both correlation matrices ([1_correlation_analysis.ipynb](1_correlation_analysis.ipynb)) and compare how well we are able to capture biological signals using each one ([1a_compare_SPELL_vs_counts_correlation.ipynb](1a_compare_SPELL_vs_counts_correlation.ipynb))