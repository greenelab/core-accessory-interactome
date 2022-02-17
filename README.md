# Identifying the interaction of core and accessory genes in *P. aeruginosa*

**Alexandra J Lee, Georgia Doing, Samuel L. Neff, Deborah A Hogan and Casey S Greene**

**April 2020**

**University of Pennsylvania**

Clinical and environmental strains of _Pseudomonas aeruginosa_ (or _P. aeruginosa_), an opportunistic pathogen that causes difficult to treat infections, have significant genomic heterogeneity including the presence of diverse accessory genes that are only present in some strains or clades.
Both core genes, which are conserved across strains, and accessory genes have been associated with traits such as biofilm formation and virulence.
Much of what we know about core and accessory gene content comes from genome analyses.
Here, we use a [newly assembled transcriptome compendium](https://www.biorxiv.org/content/10.1101/2022.01.24.477642v1) to analyze the transcriptional patterns of core and accessory gene expression in PAO1 and PA14 strains across thousands of samples from hundreds of distinct experiments.
We found that a subset of core genes was transcriptionally stable across strain PAO1 and PA14 strain types and that these genes had fewer accessory genes with correlated expression patterns than did less stable core genes.

## Directory Structure
| Folder | Description |
| --- | --- |
| [0_explore_data](0_explore_data) | This folder contains analysis notebooks to visualize the expression data to get a sense for the variation contained.|
| [1_processing](1_processing) | This folder contains analysis notebooks to determine what threshold to use to partition the gene expression data into PAO1 and PA14 compendia.|
| [2_correlation_analysis](2_correlation_analysis) | This folder contains analysis notebooks to detect gene co-expression modules starting with gene expression data, applying Pearson correlation and then clustering on this correlation matrix to obtain gene modules.|
| [3_core_core_analysis](3_core_core_analysis) | This folder contains analysis notebooks to examine the stability of core genes across strains.|
| [4_acc_acc_analysis](4_acc_acc_analysis) | This folder contains analysis notebooks to examine accessory-accessory gene modules.|
| [5_core_acc_analysis](5_core_acc_analysis) | This folder contains analysis notebooks to examine the relationship between core genes and accessory genes.|
| [6_common_genes_analysis](6_common_genes_analysis) | This folder contains analysis notebooks to compare common DEGs found in [prior work](https://github.com/greenelab/generic-expression-patterns/blob/master/pseudomonas_analysis/2_identify_generic_genes_pathways.ipynb) to core and accessory genes|
| [scripts](scripts) | This folder contains supporting functions that other notebooks in this repository will use.|
| [data](data) | This folder contains metadata used for different analyses.|


## Usage
*Operating Systems:* Mac OS, Linux (Note: bioconda libraries not available in Windows)

In order to run this simulation on your own gene expression data the following steps should be performed:

First you need to set up your local repository:
1. Download and install [github's large file tracker](https://git-lfs.github.com/). Once downloaded and installed, setup git lfs by running `git lfs install`
2. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
3. Navigate to the location where you'd like the code to live and clone the `core-accessory-interactome` repository by running the following command in the terminal:
```
git clone https://github.com/greenelab/core-accessory-interactome.git
```
Note: Git automatically detects the LFS-tracked files and clones them via http.
4. Navigate into the cloned repo by running the following command in the terminal:
```
cd core-accessory-interactome
```
5. Set up your conda environment by running the following command in the terminal:
```bash
bash install.sh
```
6. Navigate to any of the analysis directories listed in the table above to see the code for how analyses were performed. To reproduce the results and figures of the paper, run the analysis directories in order.

## Acknowledgements
We would like to thank Jake Crawford for very insightful discussions about methods and interpretation of gene correlation analyses. We would also like to thank all other members of Greene lab (Natalie Davidson, Ben Heil, Ariel Hippen, David Nicholson,  Milton Pividori,  Halie Rando, Taylor Reiter) for helpful comments and code review.
