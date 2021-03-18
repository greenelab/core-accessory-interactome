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

# # Test references created
#
# Here we want to see how well _pilA_ gene aligns against PAO1+phage reference vs PA14+phage reference. This _pilA_ gene has a low sequence identity with ????

import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from core_acc_modules import paths

# + magic_args="-s $paths.PAO1_PHAGE_REF $paths.PAO1_PHAGE_DB_DIR" language="bash"
#
# makeblastdb -in $1 -dbtype nucl -parse_seqids -out $2

# + magic_args="-s $paths.PA14_PHAGE_REF $paths.PA14_PHAGE_DB_DIR" language="bash"
#
# makeblastdb -in $1 -dbtype nucl -parse_seqids -out $2

# + magic_args="-s $paths.PILA_QUERY $paths.PAO1_PILA_BLAST_RESULT $paths.PAO1_PHAGE_DB_DIR" language="bash"
# blastn -query $1 -out $2 -db $3 -outfmt 6

# + magic_args="-s $paths.PILA_QUERY $paths.PA14_PILA_BLAST_RESULT $paths.PA14_PHAGE_DB_DIR" language="bash"
# blastn -query $1 -out $2 -db $3 -outfmt 6
# -

pao1_blast_result = pd.read_csv(paths.PAO1_PILA_BLAST_RESULT, sep="\t", header=None)
pa14_blast_result = pd.read_csv(paths.PA14_PILA_BLAST_RESULT, sep="\t", header=None)

# Add column names described above
col_names = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

# BLAST results for PAO1
pao1_blast_result.columns = col_names
print(pao1_blast_result.shape)
print(pao1_blast_result["evalue"].max())
pao1_blast_result.head()

# BLAST results for PA14
pa14_blast_result.columns = col_names
print(pa14_blast_result.shape)
print(pa14_blast_result["evalue"].max())
pa14_blast_result.head()
