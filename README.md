# Identifying the interaction of core and accessory genes in *P. aeruginosa*

**Alexandra J Lee, Georgia Doing, Deborah A Hogan and Casey S Greene**

**April 2020**

**University of Pennsylvania**


*Motivation:*  In general, it has been suggested that *P aeruginosa* adaptive traits, like virulence, depends on both core and accessory genes. As an example, one study by [Vasquez-Rifo et. al](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1890-1) found that accessory genes are statistically associated with virulence. Another study by  [Tan et. al](https://www.pnas.org/content/pnas/96/2/715.full.pdf) found that mutant gacA genes were less virulence to *C. elegans* compared to WT, suggesting a role of the global gacA regulator in virulence. Overall, adaptive traits, like virulence, appear to be multifactorial. 


*Objective:* Given that adaptive traits, like virulence, appear to be the result of multiple genes and that both core and accessory genes contribute to these traits, we want to know what the relationship is between these two sets of genes, which is currently unclear.

*Hypothesis:* Flexible genes have analogous expression with core genes than unique genes. In other words we would expect flexible and core genes to be co-expressed compared to unique and core genes. The intuition for this hypothesis was suggested by [Jiao et. al](https://www.ncbi.nlm.nih.gov/pubmed/29795552),  who found that more conserved genes are more highly connected in the co-expression network in *S. fredii* bacteria

*Expected Outcome:* We expect to have define the co-expression between core and the accessory genome.   
Knowing this information will allow us to integrate data to be able to identify possibly novel interactions instead of analyzing data in buckets. Futhermore, we can start to ask questions about the regulation pattern between core and accessory, which can further help to inform mechanisms related to *P. aeruginosa* adaptation.

## Computational Environment

All processing and analysis scripts were performed using the conda environment specified in `environment.yml`.
To build and activate this environment run:

```bash
# conda version 4.6.12
conda env create -f environment.yml

conda activate core_acc
```