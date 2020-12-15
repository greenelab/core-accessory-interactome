# Pilot: exploring core-accessory relationships within an individual experiment

**April 2020**

This pilot analysis will serve to answer the question: What is the simplest analysis that will convince us that [aim 2](../README.md) is worth pursuing?


*Question:*  What relationships between core and accessory genes can we find using small subset of samples? 

*Approach:*
1. Select an experiment from the [*P. aeruginosa* compendium](https://msystems.asm.org/content/1/1/e00025-15)
2. Map PAO1-only, PA14-only, PAO1/PA14 homolog annotations to transcriptome data
Note: We are using PAO1 and PA14 genotypes since they are well established.
3. Explore core-accessory relationships: What is the correlation between core-only genes, accessory-only genes, core-accessory genes? What explains the correlation structure (i.e. co-operonic, shared pathways)? What is known about core-accessory relationships?

*Conclusion:*
* We observed that there is a higher likelihood of accessory genes being higher correlated with other accessory genes (using PAO1 only samples). However it is difficult to interpret the results due to the limitations of using microarrays.
* With array data, we are limited to detecting only those genes that are accessory to the reference transcriptome that is used for the arrays, in this case we can detect genes that are PAO1-specific, but not PA14-specific. So if we use rna-seq we can map sequences to both PAO1 and PA14 references and use core genes = homologs between PAO1 and PA14. And accessory genes for PAO1, PA14 are all remaining genes that are not core. 
* If we look at the expression of a gene that is absent in PA14 in PA14 samples only, we can still get nonzero expression because of non-specific binding (based on measurements). There seems to be a fair amount of cross-hybridization for genes that are annotated as not having PA14 homologs (what we’re calling accessory genes) but showing signal in PA14 samples. Using RNAseq alignment parameters may fix this or at least provide clarity as to whether it is truly non-specific (random biochemical events) or whether there are genes that are similar in sequence enough to hybridized or allign (with an error rate) but not to be called homologs
* If we look at the expression of a gene that is absent in PA14 in PA14 samples only, we can still get nonzero expression even if there was NO non-specific binding, because of how the data is processed, where the variation will be attributed to the noise term in the regression. It will be important to distinguish noise from sequence similarity that leads to cross-hybridization. Reads will let us do that array signal does not. The processing and normalization is more straight forward I think and less likely to lead to 'false positives' due to normalization.
