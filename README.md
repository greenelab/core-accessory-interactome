# Identifying the interaction of core and accessory genes in *P. aeruginosa*

**Alexandra J Lee, Georgia Doing, Deborah A Hogan and Casey S Greene**

**April 2020**

**University of Pennsylvania**

_Background_:

Pseudomonas aeruginosa (or _P. aeruginosa_) is a gram negative bacteria. _P. aeruginosa_ is the most common gram-negative pathogen causing nosocomial (NO SO comial) pneumonia in the United States, and it is frequently implicated in other hospital-acquired infections like urinary tract and bloodstream infections ([ National Nosocomial Infections Surveillance System](https://academic.oup.com/cid/article/41/6/848/2022258)). In additional to its prevalence, _P. aeruginosa_ is able to develop resistance to antibiotics, which is the standard of care for these infections. Overall the predilection for _P. aeruginosa_ to cause infections in immunocompromised individuals, its extreme versatility (cause infection across different tissues) and antibiotic resistance makes _P. aeruginosa_ a major health concern.

Over the years, there have been are many studies trying to understand the genetic mechanisms that drive _P. aeruginosa_ infections. Some of these studies revealed that the Pa genome was comprised of two components: core genome and an accessory genome. Core genes are those genes that are present in all strains
While accessory genes are those that are present in at least one strain. The current paradigm is that the core genes encode essential functions shared by all strains while accessory genes encode niche-based adaptations

In general, genomic analyses of pseudomonas usually separate between core genes and accessory genes. However, in general, there is evidence that both core and accessory gene expression can contribute to phenotypes. One example is Pa virulence ([Lee et. al. 2014](https://pubmed.ncbi.nlm.nih.gov/25249263/), [Fink-Barbancon et. al. 1997](https://pubmed.ncbi.nlm.nih.gov/9302017/))

_Central question_: How is the expression of these gene groups coordinated to drive phenotypes of interest?

To answer this central question we first need to examine the relationship between genes and function, as an intermediate step before we get to phenotype. So for this project we will focus on first answering the question: **How are functions carried out in the context of core and accessory genes?**

_Hypothesis_: A common subset of core genes that will copt with accessory genes across strains

_Expected Outcome_: The composition of network modules in the context of core, accessory

_Significance_: If we know the composition of modules this can inform how we study functions. For example, we can build models to combine different gene groups and use these models to study phenotypes, which may reveal novel combinations of genes


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