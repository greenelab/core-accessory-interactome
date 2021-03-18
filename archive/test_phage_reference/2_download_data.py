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
#     display_name: Python [conda env:core_acc] *
#     language: python
#     name: conda-env-core_acc-py
# ---

# # Download data
#
# This notebook is downloading data from SRA and then using the module, `fasterq-dump`, from the SRA toolkit to get the fastq files associated with the downloaded SRA files.
#
# Note: Need to delete `sra` folder between runs; otherwise, `fastq-dump` will be called on all files in `sra` folder which can include more than your sra accessions.

# +
# %load_ext autoreload
# %autoreload 2

import os
import shutil
import pandas as pd
import numpy as np
from core_acc_modules import paths

np.random.seed(123)
# -

# ### Download SRA data
#
# Note: Need to delete `sra` folder between runs otherwise `fastq-dump` will be called on all files in `sra` folder which can include more than your sra accessions.

shutil.rmtree(paths.SRA_DIR)

# Download sra data files
# ! prefetch --option-file $paths.SRA_ACC

# ### Get FASTQ files associated with SRA downloads
#
# The fastq files store the RNA-seq results, including: sequencing and quality scores for each base call.
#
# Here is a nice blog to explain how to read fastq files: https://thesequencingcenter.com/knowledge-base/fastq-files/
#
# The fastq files gives the sequence of a read at a given location. Our goal is to map these reads to a reference genome so that we can quantify the number of reads that are at a given location, to determine the level of expression.
#
# `fasterq-dump` automatically splits paired-end data into 3 files:
# * file_1.fastq having read 1
# * file_2.fastq having read 2
# * file.fastq having unmatched reads (i.e. read doesn't have a mate pair).
# https://www.rdocumentation.org/packages/geomedb/versions/2.0.1/topics/fasterqDump

os.makedirs(paths.FASTQ_DIR, exist_ok=True)

# + magic_args="-s $paths.SRA_DIR $paths.FASTQ_DIR" language="bash"
# for FILE_PATH in $1/*;
#
# do
#     # get file name
#     sample_name=`basename ${FILE_PATH}`
#
#     # remove extension from file name
#     sample_name="${sample_name%.*}"
#
#     echo "Getting fastq for $FILE_PATH"
#     fasterq-dump $FILE_PATH -O $2/${sample_name} -f
# done
