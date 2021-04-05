# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1+dev
#   kernelspec:
#     display_name: Python [conda env:core_acc_env] *
#     language: python
#     name: conda-env-core_acc_env-py
# ---

# # Decide thresholds
#
# This notebook is examining the distribution of the mapping rates to be used for binning.

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from textwrap import fill
from core_acc_modules import paths

# Log files
pao1_logs_filename = paths.PAO1_LOGS
pa14_logs_filename = paths.PA14_LOGS

pao1_logs = pd.read_csv(pao1_logs_filename, index_col=0, header=0)
pa14_logs = pd.read_csv(pa14_logs_filename, index_col=0, header=0)

pao1_logs.head()

pa14_logs.head()

# +
# Plot distribution of mapping rates to PAO1 and PA14

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8, 4))

# Distribution plot for core genes
sns.distplot(
    pao1_logs["mapping_rate"],
    label="PAO1 mapping rate",
    color="red",
    kde=False,
    ax=axes[0],
)

sns.distplot(
    pa14_logs["mapping_rate"],
    label="PA14 mapping rate",
    color="blue",
    kde=False,
    ax=axes[1],
)

plt.suptitle(fill("Distribution of mapping rates", width=40), x=0.5, y=1, fontsize=16)
axes[0].set_title(fill("PAO1 mapping rate", width=20))
axes[1].set_title(fill("PA14 mapping rate", width=20))
axes[0].set_xlabel("")
axes[1].set_xlabel("")
fig.text(0.5, 0.01, "Mapping rate", ha="center", fontsize=14)
fig.text(0.01, 0.5, "Count", ha="center", rotation=90, fontsize=14)

# +
# Plot distribution of difference mapping rates to PAO1 and PA14

# Set up the matplotlib figure
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8, 4))

# Distribution plot for core genes
sns.distplot(
    pao1_logs["mapping_rate"] - pa14_logs["mapping_rate"],
    label="PAO1 mapping rate",
    color="red",
    kde=False,
    ax=axes[0],
)

sns.distplot(
    pa14_logs["mapping_rate"] - pao1_logs["mapping_rate"],
    label="PA14 mapping rate",
    color="blue",
    kde=False,
    ax=axes[1],
)

plt.suptitle(
    fill("Distribution of difference in mapping rates", width=40),
    x=0.5,
    y=1.1,
    fontsize=16,
)
axes[0].set_title(fill("PAO1-PA14 mapping rate", width=20))
axes[1].set_title(fill("PA14-PAO1 mapping rate", width=20))
axes[0].set_xlabel("")
axes[1].set_xlabel("")
fig.text(0.5, 0.01, "Difference in mapping rate", ha="center", fontsize=14)
fig.text(0.01, 0.5, "Count", ha="center", rotation=90, fontsize=14)
# -

# **Observations:**
# * There is fairly rough bimodal distribution, so most samples align well to PAO1 reference or the PA14 reference.
#
# **Takeaway:**
# * Based on the mapping rate distribution, looks like a mapping rate of 30% would be a good cutoff to identify those samples that have a high mapping to PAO1 or PA14
# * Since PAO1 and PA14 have a larged shared core genome, we need to add an additional criteria to filter for those that are specific to PAO1 versus PA14. This additional criteria will use the difference in mapping rates
#
# A sample will be PAO1 if:
# 1. PAO1 mapping rate >= 30%
# 2. PAO1-PA14 mapping rate > 0%
#
# Similarly for PA14 samples.
