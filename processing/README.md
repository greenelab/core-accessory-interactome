# Processing
All samples were quantified in Salmon using both a PAO1 and PA14 references.

To determine which samples are PAO1 versus PA14 we will use the median expression of accessory genes to determine if a sample is PAO1 or PA14.
A sample is considered PAO1 if the median gene expression of PA14 accessory genes is 0 and PAO1 accessory genes in > 0.
Similarlty, a sample is considered PA14 if the median gene expression of PA14 accessory genes is > 0 and PAO1 accessory genes in 0.

See "TPM of accesory genes" plot in [exploratory notebook](../explore_data/cluster_by_accessory_gene.ipynb)