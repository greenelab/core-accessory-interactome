# Identifying the interaction of core and accessory genes in *P. aeruginosa*

**Alexandra J Lee, Georgia Doing, Deborah A Hogan and Casey S Greene**

**April 2020**

**University of Pennsylvania**

_Pseudomonas aeruginosa_ (or _P. aeruginosa_) is a gram negative bacteria. _P. aeruginosa_ is the most common gram-negative pathogen causing nosocomial pneumonia in the United States, and it is frequently implicated in other hospital-acquired infections like urinary tract and bloodstream infections
 ([National Nosocomial Infections Surveillance System](https://academic.oup.com/cid/article/41/6/848/2022258)).
 In addition to its prevalence, _P. aeruginosa_ is able to develop resistance to antibiotics, which is the standard of care for these infections.
 These infections are a major concern to hospitalized patients, as they are found to be correlated with poor prognosis. Overall the predilection for _P. aeruginosa_ to cause infections in immunocompromised individuals, its extreme versatility (cause infection across different tissues) and antibiotic resistance make Pa a major health concern.

To combat these _P. aeruginosa_ infections, there have been numerous transcriptomic studies trying to understand how Pseudomonas aeruginosa genes influence traits like virulence and pathogenicity.
The hope is that these studies will help us to develop better treatment options. These contributing genes can be classified into two groups: core and accessory.
Core genes are those genes that are present in all strains.
While accessory genes are those that are present in at least one strain.

Given that different groups of genes both contribute to these traits of interest, it is important to understand how these gene groups are coordinated.

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