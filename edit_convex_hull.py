from __future__ import division
from Bio.PDB import *
import pandas as pd
from util_func import conv_array_text, extract_aa_sequence, extract_all_atoms, extract_backbone_atoms
from calc_func import convex_hull_peeling

p = PDBParser()

#  AlphaFold structure for test
structure = p.get_structure('AF-O60361-F1-model_v4',
                            r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\AF-O60361-F1-model_v4.pdb")

# Extract amino acid sequence from the first model and chain
print(extract_aa_sequence(structure))

# Extract all atomic coordinates
print(extract_all_atoms(structure))

# Extract backbone atomic coordinates
points = extract_backbone_atoms(structure)
print(points)

# Layer the structure using convex hull peeling method
convex_hull_peeling(points)

# Extract all peptide cleavage XYZ coordinates (based on time point)
data_df = pd.read_csv("pep_cleave_coordinates_10292023.csv", index_col=0)
data_df = data_df.applymap(conv_array_text)
print(data_df)

# Extract specific time-lapsed peptide cleave XYZ coordinates for a given UniProt ID (ex: O60361)
# Protein ID should match with that of the AlphaFold file
search_value = "O60361"
peptide_coord = data_df.loc[search_value]
print(peptide_coord)
