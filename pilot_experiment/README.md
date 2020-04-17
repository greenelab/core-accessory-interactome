# Pilot: exploring core-accessory relationships within an individual experiment

**April 2020**

This pilot analysis will serve to answer the question: What is the simplest analysis that will convince us that [aim 2](../README.md) is worth pursuing?


*Question:*  What relationships between core and accessory genes can we find using an individual experiment? 

*Approach:*
1. Select an experiment from the [*P. aeruginosa* compendium](https://msystems.asm.org/content/1/1/e00025-15)
2. Map PAO1-only, PA14-only, PAO1/PA14 homolog annotations to transcriptome data
Note: We are using PAO1 and PA14 genotypes since they are well established.
3. Use [WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559) to create gene-gene co-expression networks using real expression data and permuted expression data (baseline)

Validate networks by ensuring that the networks capture known operons between core and accessory genes

Explore core-accessory relationships:
1. Select pair of novel core, accessory genes that are interacting
2. Search literature for what is known about core, accessory genes. 
3. Are core and accessory genes found in the same pathways? have related functions?
4. What does the gene expression pattern look like (core vs accessory)? What is the experimental context associated with well correlated vs not well correlated samples?

*Expected Outcome:* We expect to find novel interactions between core and accessory genes that have the potential to be important for *P. aeruginosa* adaptation based on the properties found in the literature. We can also start to speculate about the role of accessory genes coordinating with core genes based on pathway association and pairwise gene expression patterns.

Based on the findings of this experiment, we might want to expand this analysis to look at:
1. What relationships change when we examine different experiments or add more experiments? Perhaps some relationships are more robust than others.
2. If we use different similarity or clustering methods, how will our relationships change? Can we detect both linear and nonlinear relationships?
3. If we use different strains, say different PAO1 strains, to define core and accessory, can we detect more subtle patterns?