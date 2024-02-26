from __future__ import division
from Bio.PDB import *
import pandas as pd
from util_func import conv_array_text,conv_array_coord, extract_aa_sequence, extract_all_atoms, extract_backbone_atoms
from calc_func import convex_hull_peeling, protein_width, average_thickness, defined_convex_hull_peeling
from graph_func import plot_xyz_coordinates_with_time_and_hulls, plot_hulls_for_time_points

p = PDBParser()

#  AlphaFold structure for test
structure = p.get_structure('Q02543',
                            r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\AF-Q02543-F1-model_v4.pdb")

# Extract amino acid sequence from the first model and chain
print(extract_aa_sequence(structure))

# Extract all atomic coordinates
print(extract_all_atoms(structure))

# Extract backbone atomic coordinates
points = extract_backbone_atoms(structure)
print(points)

# Layer the structure using convex hull peeling method
hull_layers = convex_hull_peeling(points)
#defined_layers = defined_convex_hull_peeling(points)

# Extract all peptide cleavage XYZ coordinates (based on time point)
data_df = pd.read_csv("pep_cleave_coordinates_10292023.csv", index_col=0)
data_df = data_df.applymap(conv_array_text)

# Extract specific time-lapsed peptide cleave XYZ coordinates for a given UniProt ID (ex: O60361)
# Protein ID should match with that of the AlphaFold file
search_value = "Q02543"
peptide_coord = data_df.loc[search_value]
print(peptide_coord)

conv_peptide_coord = conv_array_coord(peptide_coord)
print(conv_peptide_coord)

# Check for breakdown pattern by indicating detected layers at each time point
plot_hulls_for_time_points(conv_peptide_coord, hull_layers)
#plot_hulls_for_time_points(conv_peptide_coord, defined_layers)


# Plot proteolytic data onto layered structure
plot_xyz_coordinates_with_time_and_hulls(conv_peptide_coord, hull_layers)
#plot_xyz_coordinates_with_time_and_hulls(conv_peptide_coord, defined_layers)


# Calculate the overall max width of the protein and average thickness of each layer (in Angstroms)
print(protein_width(points))
print(average_thickness(hull_layers))

def detect_aa_residues_in_hulls(structure, hull_layers):
    amino_acid_sequence = extract_aa_sequence(structure)
    lysine_layers = []
    arginine_layers = []

    for i, hull_points in enumerate(hull_layers):
        lysine_count = 0
        arginine_count = 0

        for residue in amino_acid_sequence:
            if residue == 'K':
                lysine_count += 1
            elif residue == 'R':
                arginine_count += 1

        lysine_layers.append(lysine_count)
        arginine_layers.append(arginine_count)

    return lysine_layers, arginine_layers
residues = detect_aa_residues_in_hulls(points, hull_layers)
print(residues)
