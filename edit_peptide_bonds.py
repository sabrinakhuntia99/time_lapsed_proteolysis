import numpy as np
import matplotlib.pyplot as plt
from util_func import extract_backbone_atoms
from Bio.PDB import *

p = PDBParser()

#  AlphaFold structure for test
structure = p.get_structure('O60361',
                            r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\AF-O60361-F1-model_v4.pdb")

def calculate_peptide_bond_rotation(structure):
    rotations = []
    for model in structure:
        for chain in model:
            prev_residue = None
            for residue in chain:
                if prev_residue is not None:
                    for atom1, atom2, atom3, atom4 in zip(prev_residue["C"], residue["N"], residue["Cα"], residue["N"]):
                        dihedral_angle = calcula           te_dihedral_angle(atom1.get_coord(), atom2.get_coord(), atom3.get_coord(), atom4.get_coord())
                        rotations.append(dihedral_angle)
                prev_residue = residue
    return rotations

def calculate_dihedral_angle(p1, p2, p3, p4):
    b1 = -1.0 * (p2 - p1)
    b2 = p3 - p2
    b3 = p4 - p3

    b2 /= np.linalg.norm(b2)

    v = b1 - np.dot(b1, b2) * b2
    w = b3 - np.dot(b3, b2) * b2

    x = np.dot(v, w)
    y = np.dot(np.cross(b2, v), w)

    return np.degrees(np.arctan2(y, x))

def assign_flexibility(rotations):
    flexibility_score = np.std(rotations)
    return flexibility_score
def plot_backbone_with_flexibility(structure, rotations):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    backbone_atoms = extract_backbone_atoms(structure)
    min_flexibility = min(rotations)
    max_flexibility = max(rotations)
    scaled_flexibility = [(flex - min_flexibility) / (max_flexibility - min_flexibility) for flex in rotations]

    ax.scatter(backbone_atoms[:, 0], backbone_atoms[:, 1], backbone_atoms[:, 2], color='black', label='Backbone Atoms')

    for i in range(len(backbone_atoms) - 4):
        color_intensity = scaled_flexibility[i]
        ax.plot(
            [backbone_atoms[i, 0], backbone_atoms[i + 1, 0]],
            [backbone_atoms[i, 1], backbone_atoms[i + 1, 1]],
            [backbone_atoms[i, 2], backbone_atoms[i + 1, 2]],
            color=(color_intensity, 0, 0),
            linewidth=2
        )

# Calculate peptide bond rotations
rotations = calculate_peptide_bond_rotation(structure)

# Assign flexibility
flexibility_score = assign_flexibility(rotations)
print("Flexibility Score:", flexibility_score)

# Plot backbone atoms with flexibility-colored peptide bonds
plot_backbone_with_flexibility(structure, rotations)