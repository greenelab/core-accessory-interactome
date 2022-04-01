# Processing
The Raw data was quantified in Salmon using both PAO1 and PA14 references.
For more information on the raw data is see this [Doing et al.](https://www.biorxiv.org/content/10.1101/2022.01.24.477642v1) with source code [here](https://github.com/hoganlab-dartmouth/pa-seq-compendia).
The datasets containing all strains aligned against the PAO1 reference and PA14 reference is [here](https://osf.io/vz42h/).

To determine which samples were more PAO1- versus PA14-like, we will use the median expression of accessory genes.
In our exploratory analysis we found that strains labeled as PAO1 based on SRA annotations had high PAO1 accessory gene expression.
Whereas strains labeled as PA14 by SRA had high PA14 accessory gene expression.

See plot below where the median expression of PAO1 genes (PAO1 accessory genes) on the x-axis and the median expression of PA14-only genes (PA14 accessory genes) on the y-axis.
Each point is a strain/sample.
![all_samples](https://github.com/greenelab/core-accessory-interactome/blob/master/0_explore_data/Expression_accessory_genes_all_samples.svg)

A strain is considered PAO1 if the median gene expression of PA14 accessory genes is 0 and PAO1 accessory genes in > 5.
Similarlty, a strain is considered PA14 if the median gene expression of PA14 accessory genes is > 5 and PAO1 accessory genes in 0.

A threshold of 25 MR normalized estimated counts is used based on our analysis in [1_decide_threshold.ipynb](1_decide_threshold.ipynb). The goal of this notebook was to  define a threshold to determine if a strain if PAO1 or not (likewise, if a strain is PA14 or not). We used known labels from SRA to do this. Specifically, we examined the distribution of PAO1 strain vs non-PAO1 strain (see histogram plots). We define the threshold to be one that separated between the two distributions. We use this threshold in [2_create_compendia.ipynb](2_create_compendia.ipynb) to partition gene expression data into PAO1 and PA14 compendia because we found that using a threshold of 0 MR normalized estimated counts included some other SRA-labeled strains.

Using a threshold of 0 MR normalized estimated counts, within the PAO1 binned compendium there are strains that SRA labeled as PAK or Clinical.
![pao1_compendium_0thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/MR_median_acc_expression_pao1_compendium_0threshold.svg)


Similarly for the PA14 compendium.

![pa14_compendium_0thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/MR_median_acc_expression_pa14_compendium_0threshold.svg)

Looking at the distribution of the median accessory gene expression for these non-PAO1 SRA labeled strains (i.e. PAK, Clinical) their expression is very low (brown) compared to all other PAO1 labeled strains (light purple/pink).
A similar trend is seen comparing the non-PA14 labeled strains (brown) vs PA14 strains (dark purple)

![pao1_dist_0thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/dist_median_acc_expression_pao1_compendium_0threshold.svg)

![pa14_dist_0thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/dist_median_acc_expression_pa14_compendium_0threshold.svg)

Note: To generate the above figures, set the params for `same_threshold = 0` and `opp_threshold = 25` in [2_create_compendia](2_create_compendia.ipynb) noetbook. You will also need to comment out the assertion statement: `assert len(shared_pao1_pa14_binned_ids) == 0`. Then run [3_validate_compendia](3_validate_compendia.ipynb).

Using a threshold of 25 we get the following plots that correspond to our final compendia that we will use in our analysis.
As a check, our PAO1 compendium contains 890 strains and the PA14 compendium contains 505 strains.
These numbers are close to the numbers that SRA annotates as PAO1 and PA14.

![pao1_compendium_25thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/MR_median_acc_expression_pao1_compendium_25threshold.svg)

![pa14_compendium_25thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/MR_median_acc_expression_pa14_compendium_25threshold.svg)

![pao1_dist_25thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/dist_median_acc_expression_pao1_compendium_25threshold.svg)

![pa14_dist_25thresdhold](https://github.com/greenelab/core-accessory-interactome/blob/master/1_processing/dist_median_acc_expression_pa14_compendium_25threshold.svg)
