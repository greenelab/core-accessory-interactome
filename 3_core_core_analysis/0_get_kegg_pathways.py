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

# +
from Bio.KEGG import REST

pao1_pathways = REST.kegg_list("pathway", "pae").read()

pao1_pathways.split("\n")[0:5]

# +
# This python package or the R getGeneSet both return only the KEGG pathways but not the modules
# Not sure how to get modules but don't worry about this right now
# Python we will need to parse each line
# OR R need to figure out how to map nonsense ids to gene ids
# Need to make sure to have pathway name and id too
# https://biopython-tutorial.readthedocs.io/en/latest/notebooks/18%20-%20KEGG.html

# +
# Clean up pathway names
pao1_pathways_clean = []
for line in pao1_pathways.rstrip().split("\n"):
    entry, description = line.split("\t")
    pao1_pathways_clean.append(entry)

print(len(pao1_pathways_clean))
pao1_pathways_clean

# +
# List of all pae pathways in KEGG
# as a check there are 123
# https://www.kegg.jp/kegg-bin/search_pathway_text?map=pae&keyword=&mode=1&viewImage=true

# +
# Example of a KEGG pathway page that is being parsed
# https://www.genome.jp/dbget-bin/www_bget?pathway:pae00020
# For a given pathway there are multiple modules and multiple gene ids

# +
# Get the genes for pathways and add them to a list
pathway_id_to_genes = {}
pathway_id_to_pathway_name = {}

for pathway in pao1_pathways_clean:
    print(pathway)
    pathway_file = REST.kegg_get(pathway).read()  # query and read each pathway

    # iterate through each KEGG pathway file, keeping track of which section
    # of the file we're in, only read the gene in each pathway

    # Initiate empty list to store genes
    pathway_gene_ids = []
    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        print(line)
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section

        if current_section == "NAME":
            pathway_id_to_pathway_name[pathway] = line

        if current_section == "GENE":
            print("line in gene", line[12:])
            if line.count(";") == 1:
                gene_identifiers, gene_description = line[12:].split("; ")
                gene_id, gene_symbol = gene_identifiers.split()
            elif line.count(";") == 0:
                gene_identifiers = line[12:].split(" ")
                gene_id = gene_identifiers[0]
            else:
                # in the case of 2 semicolons
                gene_identifiers, gene_description, protein_name = line[12:].split("; ")
                gene_id, gene_symbol = gene_identifiers.split()

            if gene_id not in pathway_gene_ids:
                pathway_gene_ids.append(gene_id)

    pathway_id_to_genes[pathway] = pathway_gene_ids
# -

pathway_id_to_genes["path:pae01110"]

# +
# import pandas as pd
# pathway_gene_df = pd.DataFrame.from_dict(pathway_id_to_genes)
# Need a two column pandas df with pathway - list

# +
# "PA2624" in pathway_id_to_genes.values()

# +
# import requests as r
# r.requests.get("http://rest.kegg.jp/list/pathway/pae")
