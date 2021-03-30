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
from core_acc_modules import paths

# Log files
pao1_logs_filename = paths.PAO1_LOGS
pa14_logs_filename = paths.PA14_LOGS

pao1_logs = pd.read_csv(pao1_logs_filename, index_col=0, header=0)
pa14_logs = pd.read_csv(pa14_logs_filename, index_col=0, header=0)

pao1_logs.head()

pa14_logs.head()

sns.distplot(pao1_logs["mapping_rate"], kde=False)
sns.distplot(pa14_logs["mapping_rate"], kde=False)

# **Observations:**
# * There is fairly rough bimodal distribution, as expected, as most samples should align well to PAO1 reference or the PA14 reference.
#
# **Takeaway:**
# * Based on the distribution, looks like a mapping rate of 25% would be a good cutoff to use for both compendia
