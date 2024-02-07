from __future__ import division
from Bio.PDB import *
import pandas as pd
from util_func import conv_array_text,conv_array_coord, extract_aa_sequence, extract_all_atoms, extract_backbone_atoms
from calc_func import convex_hull_peeling, calculate_average_depth_of_hulls
from graph_func import plot_xyz_coordinates_with_time_and_hulls, plot_hulls_for_time_points

p = PDBParser()

#  AlphaFold structure for test
structure = p.get_structure('O60361',
                            r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\AF-O60361-F1-model_v4.pdb")

# Extract amino acid sequence from the first model and chain
print(extract_aa_sequence(structure))

# Extract all atomic coordinates
print(extract_all_atoms(structure))

# Extract backbone atomic coordinates
points = extract_backbone_atoms(structure)
print(points)

# Layer the structure using convex hull peeling method
hull_layers = convex_hull_peeling(points)
average_depth = calculate_average_depth_of_hulls(hull_layers)
print("Average Depth of Convex Hull Layers:", average_depth)

# Extract all peptide cleavage XYZ coordinates (based on time point)
data_df = pd.read_csv("pep_cleave_coordinates_10292023.csv", index_col=0)
data_df = data_df.applymap(conv_array_text)

# Extract specific time-lapsed peptide cleave XYZ coordinates for a given UniProt ID (ex: O60361)
# Protein ID should match with that of the AlphaFold file
search_value = "O60361"
peptide_coord = data_df.loc[search_value]
print(peptide_coord)

conv_peptide_coord = conv_array_coord(peptide_coord)
print(conv_peptide_coord)

# Plot proteolytic data onto layered structure
plot_xyz_coordinates_with_time_and_hulls(conv_peptide_coord, hull_layers)

# Check for breakdown pattern by indicating detected layers at each time point
plot_hulls_for_time_points(conv_peptide_coord, hull_layers)
