from __future__ import division
from Bio.PDB import *
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.spatial import distance
from scipy.spatial import ConvexHull
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
from util_func import conv_array_text, extract_backbone_atoms
from calc_func import convex_hull_peeling
from tqdm import tqdm

parse = PDBParser()

# Extract all peptide cleavage XYZ coordinates (based on time point)
data_df = pd.read_csv("pep_cleave_coordinates_10292023.csv", index_col=0)
data_df = data_df.applymap(conv_array_text)
data_df = data_df[0:5]

# Atomic masses in Dalton (g/mol)
atomic_masses = {
    'H': 1.008,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'P': 30.974,
    'S': 32.06,
    'SE': 78.96,
}

# Preallocate arrays for storing disorder_properties and CPI
num_proteins = len(data_df.index)
disorder_properties = {
    'Density': np.zeros(num_proteins),
    'Radius_of_Gyration': np.zeros(num_proteins),
    'Surface_Area_to_Volume_Ratio': np.zeros(num_proteins),
    'Sphericity': np.zeros(num_proteins),
    'CPI': np.zeros(num_proteins)
}

def calculate_density(structure):
    atoms = list(structure.get_atoms())
    total_mass = sum(atomic_masses.get(atom.element, 0) for atom in atoms)

    atoms_coords = [atom.coord for atom in atoms]
    hull = ConvexHull(atoms_coords)
    volume = hull.volume  # in cubic Angstroms

    # Convert volume to cubic centimeters (1 Å³ = 1e-24 cm³)
    volume_cm3 = volume * 1e-24

    # Convert mass from Daltons to grams (1 Da = 1.66054e-24 g)
    mass_g = total_mass * 1.66054e-24

    # Calculate density in g/cm³
    density = mass_g / volume_cm3
    return density

def calculate_radius_of_gyration(structure):
    atoms_coords = np.array([atom.coord for atom in structure.get_atoms()])
    center_of_mass = np.mean(atoms_coords, axis=0)
    distances = distance.cdist(atoms_coords, [center_of_mass])
    radius_of_gyration = np.sqrt(np.mean(distances**2))
    return radius_of_gyration

def calculate_surface_area_to_volume_ratio(structure):
    atoms_coords = [atom.coord for atom in structure.get_atoms()]
    hull = ConvexHull(atoms_coords)
    surface_area = hull.area
    volume = hull.volume
    return surface_area / volume

def calculate_sphericity(structure):
    atoms_coords = [atom.coord for atom in structure.get_atoms()]
    hull = ConvexHull(atoms_coords)
    volume = hull.volume
    surface_area = hull.area
    equivalent_radius = (3 * volume / (4 * np.pi))**(1/3)
    return (np.pi**(1/3)) * ((6 * equivalent_radius)**(2/3)) / surface_area

def get_pdb_file_paths(folder_path):
    pdb_paths = {}
    pattern = re.compile(r"AF-(\w+)-F\d+-model_v4.pdb")

    for subdir, _, files in tqdm(os.walk(folder_path)):
        match = pattern.match(os.path.basename(subdir))
        if match:
            uniprot_id = match.group(1)
            pdb_files = [f for f in files if f.endswith('.pdb')]
            if pdb_files:
                pdb_paths[uniprot_id] = os.path.join(subdir, pdb_files[0])

    return pdb_paths

# Get dictionary of PDB file paths
pdb_paths_dict = get_pdb_file_paths(r"C:\Users\Sabrina\PycharmProjects\intrinsic_disorder\proteome_human")

# Create a list to store hull layer data
hull_layers_data = []

# Iterate through each UniProt ID in the data file
for uniprot_id in data_df.index:
    try:
        # Get corresponding PDB file path
        pdb_file_path = pdb_paths_dict.get(uniprot_id)

        # Check if the PDB file exists
        if not os.path.isfile(pdb_file_path):
            print(f"PDB file not found for UniProt ID {uniprot_id}. Skipping...")
            continue

        # Parse PDB file
        structure = parse.get_structure(uniprot_id, pdb_file_path)

        # Extract backbone atomic coordinates
        points = extract_backbone_atoms(structure)

        # Layer the structure using convex hull peeling method
        hull_layers = convex_hull_peeling(points)

        # Iterate over each time point header in the data file
        for timepoint_header in data_df.columns[1:]:
            timepoint = timepoint_header.split('_')[-1]

            try:
                for coord in data_df.loc[uniprot_id, timepoint_header]:
                    for i, hull_layer in enumerate(hull_layers):
                        if coord in hull_layer:
                            hull_layers_data.append({
                                'uniprot_id': uniprot_id,
                                'timepoint_header': timepoint_header,
                                'layer': f"Hull_Layer_{i + 1}",
                                'value': True
                            })

                # Mark if coordinate not found in any hull layer
                for i in range(1, 11):  # Assuming 10 hull layers
                    if f"Hull_Layer_{i}" not in data_df.columns:
                        hull_layers_data.append({
                            'uniprot_id': uniprot_id,
                            'timepoint_header': timepoint_header,
                            'layer': f"Not_In_Hull_Layer_{i}",
                            'value': True
                        })

            except Exception as e:
                print(f"An error occurred for UniProt ID {uniprot_id}: {e}")
                continue

    except Exception as e:
        continue

# Convert the list of dictionaries to a DataFrame
hull_layers_df = pd.concat([pd.DataFrame.from_records(hull_layers_data).set_index(['uniprot_id', 'timepoint_header', 'layer'])['value'].unstack(level=-1)], axis=1)
print(hull_layers_df)

# Calculate average hull layer number for each protein
average_hull_layer = hull_layers_df.mean(axis=1)
print(average_hull_layer)

# Initialize CPI column in data_df
data_df['CPI'] = np.nan

# Output the index and columns of data_df for debugging
print("Index of data_df:", data_df.index)
print("Columns of data_df:", data_df.columns)

# Assign CPI based on average hull layer number
for uniprot_id, avg_layer in average_hull_layer.items():
    if avg_layer >= 8:
        data_df.loc[uniprot_id, 'CPI'] = -1  # Outer layers only
    elif avg_layer >= 4:
        data_df.loc[uniprot_id, 'CPI'] = 0   # Middle layers
    else:
        data_df.loc[uniprot_id, 'CPI'] = 1   # Inner layers

# Create a DataFrame to store the presence of proteins in each hull layer for each timepoint
presence_df = pd.DataFrame(index=data_df.index, columns=data_df.columns[1:])

# Iterate over each time point header in the data file
for timepoint_header in data_df.columns[1:]:
    # Ensure we're only considering the timepoint headers, not the 'Hull_Layer' level
    if isinstance(timepoint_header, tuple):
        timepoint = timepoint_header[0].split('_')[-1]

        # Iterate over each protein
        for uniprot_id in data_df.index:
            try:
                presence = [False] * 10  # Initialize all hull layers as not detected

                # Check if the protein is detected in any hull layer for the current timepoint
                for hull_layer in range(1, 11):  # Hull layers are from 1 to 10
                    if hull_layers_df.loc[uniprot_id, (timepoint_header, f"Hull_Layer_{hull_layer}")]:
                        presence[hull_layer - 1] = True  # Mark hull layer as detected

                # Store the presence information in the presence DataFrame
                presence_df.loc[uniprot_id, timepoint_header] = presence

            except Exception as e:
                print(f"An error occurred for UniProt ID {uniprot_id}: {e}")
                continue

# Convert presence information to numeric (0 for False, 1 for True)
presence_df = presence_df.applymap(lambda x: [1 if y else 0 for y in x])

# Loop through each protein to calculate disorder_properties and CPI
for idx, uniprot_id in enumerate(data_df.index):
    try:
        # Get corresponding PDB file path
        pdb_file_path = pdb_paths_dict.get(uniprot_id)

        # Check if the PDB file exists
        if not os.path.isfile(pdb_file_path):
            print(f"PDB file not found for UniProt ID {uniprot_id}. Skipping...")
            continue

        # Parse PDB file
        structure = parse.get_structure(uniprot_id, pdb_file_path)

        # Calculate disorder_properties
        disorder_properties['Density'][idx] = calculate_density(structure)
        disorder_properties['Radius_of_Gyration'][idx] = calculate_radius_of_gyration(structure)
        disorder_properties['Surface_Area_to_Volume_Ratio'][idx] = calculate_surface_area_to_volume_ratio(structure)
        disorder_properties['Sphericity'][idx] = calculate_sphericity(structure)
        disorder_properties['CPI'][idx] = data_df.loc[uniprot_id, 'CPI']

    except Exception as e:
        print(f"An error occurred for UniProt ID {uniprot_id}: {e}")
        continue

# Convert disorder_properties dictionary to DataFrame
disorder_properties_df = pd.DataFrame(disorder_properties)

# Calculate correlation coefficients
correlation_coefficients = {}
for property_name in ['Density', 'Radius_of_Gyration', 'Surface_Area_to_Volume_Ratio', 'Sphericity']:
    correlation_coefficient, _ = pearsonr(disorder_properties_df[property_name], disorder_properties_df['CPI'])
    correlation_coefficients[property_name] = correlation_coefficient

# Print correlation coefficients
print("Correlation Coefficients:")
for property_name, coefficient in correlation_coefficients.items():
    print(f"{property_name}: {coefficient}")

# Plot the presence heatmap
plt.figure(figsize=(14, 8))
sns.heatmap(presence_df.transpose(), cmap="YlGnBu", cbar=False, annot=False)
plt.xlabel('Protein')
plt.ylabel('Timepoints')
plt.title('Hull Layer Presence of Proteins Over Time')
plt.show()
