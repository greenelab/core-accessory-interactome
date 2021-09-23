"""
Author: Alexandra Lee
Date Created: 10 August 2021

Scripts to analyze gene co-expression modules
"""
import pandas as pd
import numpy as np
from itertools import product


def get_intra_module_dist(annot_df, pa_prefix):
	"""
	For genes in the same module, calculate the pairwise distance from each other
	calculate the median pairwise distance to represent how spread the module is
	across the genome.

	Arguments
	----------
	annot_df (pandas dataframe): dataframe gene id x features

	pa_prefix (string): "PA" for PAO1 strains or "PA14_" for PA14 strains

	Returns
	--------
	Dataframe that is gene id x [median pairwise distance, min pairwise distance, max pairwise distance]
	"""
	rows = []
	for grp_name, grp_df in annot_df.groupby("module id"):
		# print("module", grp_name)

		# Trim off "PA" and convert number to integer
		ids = grp_df.index

		# Convert trailing id numbers to floats
		num_ids = [float(_id.split(pa_prefix)[1]) for _id in ids]

		abs_dist = []
		# So `num_ids` is a list of the the gene ids, which we're
		# using as our genomic location.
		# Say our list is ["PA0001", "PA0004", "PA0010"], then product
		#  will return an iterator where at the first iteration gene1 = "PA0001"
		# and gene2="PA0001" so there is a distance of 0 here.
		# But in the next iteration we have gene1 = "PA0001" and gene2 = "PA0004"
		# so the distance is 3 here.
		for gene1, gene2 in product(num_ids, num_ids):
			if gene1 != gene2:
				dist = abs(gene1 - gene2)
				# print(gene1, gene2, dist)
				abs_dist.append(dist)

		median_module_dist = np.median(abs_dist)
		mean_module_dist = np.mean(abs_dist)
		min_dist = np.min(abs_dist)
		max_dist = np.max(abs_dist)
		range_module = max_dist - min_dist

		for _id in ids:
			rows.append(
				{
					"gene id": _id,
					"mean pairwise dist": mean_module_dist,
					"median pairwise dist": median_module_dist,
					"min pairwise dist": min_dist,
					"max pairwise dist": max_dist,
					"range pairwise dist": range_module
				}
			)

	module_dist = pd.DataFrame(rows)
	module_dist = module_dist.set_index("gene id")

	return module_dist
