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

# Extract all atoms
print(extract_all_atoms(structure))

# Extract backbone atoms
points = extract_backbone_atoms(structure)
print(points)

# Extract peptide cleavage XYZ coordinates
data_df = pd.read_csv("pep_cleave_coordinates_10292023.csv", index_col=0)
data_df = data_df.applymap(conv_array_text)
print(data_df)

# Peel convex hulls
convex_hull_peeling(points)



