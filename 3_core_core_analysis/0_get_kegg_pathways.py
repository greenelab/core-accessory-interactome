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

from Bio.KEGG import REST
import pandas as pd

pao1_pathways = REST.kegg_list("pathway", "pae").read()
pa14_pathways = REST.kegg_list("pathway", "pau").read()

# Example
pao1_pathways.split("\n")[0:5]

# ### Clean up pathway names

# +
pao1_pathways_clean = []
for line in pao1_pathways.rstrip().split("\n"):
    entry, description = line.split("\t")
    pao1_pathways_clean.append(entry)

print(len(pao1_pathways_clean))
pao1_pathways_clean

# +
pa14_pathways_clean = []
for line in pa14_pathways.rstrip().split("\n"):
    entry, description = line.split("\t")
    pa14_pathways_clean.append(entry)

print(len(pa14_pathways_clean))
pa14_pathways_clean


# +
# List of all pae pathways in KEGG
# as a check there are 123
# https://www.kegg.jp/kegg-bin/search_pathway_text?map=pae&keyword=&mode=1&viewImage=true

# +
# Example of a KEGG pathway page that is being parsed
# https://www.genome.jp/dbget-bin/www_bget?pathway:pae00020
# For a given pathway there are multiple modules and multiple gene ids
# -

# ### Get associated genes

# +
# Get the genes for pathways and add them to a list


def get_geneset(pathway_id_list):
    pathway_id_to_genes = {}
    pathway_id_to_pathway_name = {}
    pathway_id_to_gene_len = {}

    for pathway in pathway_id_list:
        pathway_file = REST.kegg_get(pathway).read()  # query and read each pathway

        # Initiate empty list to store genes
        pathway_gene_ids = []

        # iterate through each KEGG pathway file, keeping track of which section
        # of the file we're in, only read the gene in each pathway
        current_section = None
        for line in pathway_file.rstrip().split("\n"):
            section = line[:12].strip()  # section names are within 12 columns
            if not section == "":
                current_section = section

            if current_section == "NAME":
                if line.count(" - ") == 1:
                    pathway_name, strain_name = line[12:].split(" - ")
                    pathway_id_to_pathway_name[pathway] = pathway_name
                else:
                    pathway_name, strain_name_1, stain_name_2 = line[12:].split(" - ")
                    pathway_id_to_pathway_name[pathway] = pathway_name

            if current_section == "GENE":
                if line.count(";") == 1:
                    gene_identifiers, gene_description = line[12:].split("; ")
                    gene_id, gene_symbol = gene_identifiers.split()
                elif line.count(";") == 0:
                    gene_identifiers = line[12:].split(" ")
                    gene_id = gene_identifiers[0]
                else:
                    # in the case of 2 semicolons
                    gene_identifiers, gene_description, protein_name = line[12:].split(
                        "; "
                    )
                    gene_id, gene_symbol = gene_identifiers.split()

                if gene_id not in pathway_gene_ids:
                    pathway_gene_ids.append(gene_id)

        pathway_id_to_genes[pathway] = pathway_gene_ids
        pathway_id_to_gene_len[pathway] = len(pathway_gene_ids)

    return (pathway_id_to_genes, pathway_id_to_pathway_name, pathway_id_to_gene_len)


# +
# Check that there are some pathways with no genes mapped
# pathway_id_to_genes["path:pae01110"]
# -

pao1_genesets_dict, pao1_pathway_names_dict, pao1_geneset_len_dict = get_geneset(
    pao1_pathways_clean
)

pa14_genesets_dict, pa14_pathway_names_dict, pa14_geneset_len_dict = get_geneset(
    pa14_pathways_clean
)

# ### Format
#
# Make gmt pathway file with index = pathway id, name and column with list of associated gene ids

assert len(pao1_pathway_names_dict) == len(pao1_genesets_dict)
assert len(pa14_pathway_names_dict) == len(pa14_genesets_dict)

# +
pao1_pathway_name_df = pd.DataFrame(
    pao1_pathway_names_dict.items(), columns=["pathway_id", "pathway_name"]
).set_index("pathway_id")
pa14_pathway_name_df = pd.DataFrame(
    pa14_pathway_names_dict.items(), columns=["pathway_id", "pathway_name"]
).set_index("pathway_id")

pao1_pathway_name_df.head()
# -

pa14_pathway_name_df.head()

pao1_pathway_gene_df = pd.DataFrame(
    pao1_genesets_dict.items(), columns=["pathway_id", "gene_ids"]
).set_index("pathway_id")
pa14_pathway_gene_df = pd.DataFrame(
    pa14_genesets_dict.items(), columns=["pathway_id", "gene_ids"]
).set_index("pathway_id")
pao1_pathway_gene_df.head()

pa14_pathway_gene_df.head()

pao1_pathway_gene_len_df = pd.DataFrame(
    pao1_geneset_len_dict.items(), columns=["pathway_id", "num_genes"]
).set_index("pathway_id")
pa14_pathway_gene_len_df = pd.DataFrame(
    pa14_geneset_len_dict.items(), columns=["pathway_id", "num_genes"]
).set_index("pathway_id")
pao1_pathway_gene_len_df.head()

# +
# Merge dfs
pao1_tmp = pao1_pathway_name_df.merge(
    pao1_pathway_gene_len_df, left_index=True, right_index=True
)
pao1_pathway_annot_df = pao1_tmp.merge(
    pao1_pathway_gene_df, left_index=True, right_index=True
)

pa14_tmp = pa14_pathway_name_df.merge(
    pa14_pathway_gene_len_df, left_index=True, right_index=True
)
pa14_pathway_annot_df = pa14_tmp.merge(
    pa14_pathway_gene_df, left_index=True, right_index=True
)

print(pao1_pathway_annot_df.shape)
pao1_pathway_annot_df.head()

# +
# Manually check that KEGG downloaded data matches what is on the KEGG website
# and that this new data is updated from what is in https://raw.githubusercontent.com/greenelab/adage/7a4eda39d360b224268921dc1f2c14b32788ab16/Node_interpretation/pseudomonas_KEGG_terms.txt
# Verified pae00020 has 28 genes in our downloaded and online versions, while ADAGE has 34 genes
# https://www.genome.jp/dbget-bin/www_bget?pathway:pae00020
# Verified that pae00071 has 34 genes in our downloaded and online versions, while ADAGE has 32 genes
# https://www.genome.jp/entry/pathway+pae00071
# pao1_pathway_annot_df.loc["path:pae00071"]
# Verified pae00072 is not found in downloaded and online versions
# -

print(pa14_pathway_annot_df.shape)
pa14_pathway_annot_df.head()

# Join pathway id and pathway name columns
pao1_pathway_annot_df["pathway_id_name"] = (
    pao1_pathway_annot_df.index + " : " + pao1_pathway_annot_df["pathway_name"]
)
pa14_pathway_annot_df["pathway_id_name"] = (
    pa14_pathway_annot_df.index + " : " + pa14_pathway_annot_df["pathway_name"]
)

# Set index
pao1_pathway_annot_df = pao1_pathway_annot_df.set_index("pathway_id_name")
pa14_pathway_annot_df = pa14_pathway_annot_df.set_index("pathway_id_name")

pao1_pathway_annot_df.head()

pa14_pathway_annot_df.head()

# Save
pao1_pathway_annot_df.to_csv("pao1_kegg_annot.tsv", sep="\t")
pa14_pathway_annot_df.to_csv("pa14_kegg_annot.tsv", sep="\t")
