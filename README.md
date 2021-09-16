# Identifying the interaction of core and accessory genes in *P. aeruginosa*

**Alexandra J Lee, Georgia Doing, Deborah A Hogan and Casey S Greene**

**April 2020**

**University of Pennsylvania**

_Pseudomonas aeruginosa_ (or _P. aeruginosa_) is a gram negative bacteria. _P. aeruginosa_ is the most common gram-negative pathogen causing nosocomial pneumonia in the United States, and it is frequently implicated in other hospital-acquired infections like urinary tract and bloodstream infections
 ([ National Nosocomial Infections Surveillance System](https://academic.oup.com/cid/article/41/6/848/2022258)).
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
| [acc_acc_analysis](acc_acc_analysis) | This folder contains analysis notebooks to examine accessory-accessory gene modules.|
| [common_genes](common_genes) | This folder contains analysis notebooks to compare common DEGs found in [prior work](https://github.com/greenelab/generic-expression-patterns/blob/master/pseudomonas_analysis/2_identify_generic_genes_pathways.ipynb) to core and accessory genes|
| [core_acc_analysis](core_acc_analysis) | This folder contains analysis notebooks to examine the relationship between core genes and accessory genes.|
| [core_acc_modules](core_acc_modules) | This folder contains supporting functions that other notebooks in this repository will use.|
| [core_core_analysis](core_core_analysis) | This folder contains analysis notebooks to examine the stability of core genes across strains.|
| [correlation_modules](correlation_modules) | This folder contains analysis notebooks to detect gene co-expression modules starting with gene expression data, applying Pearson correlation and then clustering on this correlation matrix to obtain gene modules.|
| [data](data) | This folder contains metadata used for different analyses.|
| [explore_data](explore_data) | This folder contains analysis notebooks to visualize the expression data to get a sense for the variation contained.|
| [processing](processing) | This folder contains analysis notebooks to determine what threshold to use to partition the gene expression data into PAO1 and PA14 compendia.|


## Usage
If you want to run notebooks from core-accessory-interactome, you need to follow these steps to create a conda environment and download the necessary data.

1. In the terminal execute the following steps to setup the conda environment:
```bash
# conda version 4.6.12
conda env create -f environment.yml
conda activate core_acc
pip install -e .
```

2. Update paths in `core_acc_modules/paths.py`

3. Download reference transcriptomes from Pseudomonas.com and store in appropriate path defined in `core_acc_modules/paths.py`

4. In the terminal launch jupyter notebook within your environment:
```bash
jupyter notebook
```

**WARNING: This section is a work in progress and will be updated as the project develops**