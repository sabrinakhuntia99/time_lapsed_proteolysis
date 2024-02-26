from __future__ import division
from scipy.spatial import ConvexHull
from Bio.PDB import *
import numpy as np
import plotly.graph_objs as go

p=PDBParser()
structure=p.get_structure('pdb_test', r"C:\Users\Sabrina\PycharmProjects\convexhull\venv\pdb_test.pdb")

dna_points = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "P":
                    dna_point = atom.get_coord()
                    dna_points.append(dna_point)
dna_points = np.array(dna_points)

prot_points = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "N":
                    prot_point = atom.get_coord()
                    prot_points.append(prot_point)
prot_points = np.array(prot_points)

points = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                point = atom.get_coord()
                points.append(point)
points = np.array(points)

# calculate centroid

# centroid = np.mean(points, axis=0)
# print(centroid)

if __name__ == "__main__":
    # complex hull
    hull = ConvexHull(points)
    vertices = np.copy(hull.vertices)
    hull_points = np.copy(points[vertices])

    # DNA hull

    dna_hull = ConvexHull(dna_points)
    dna_vertices = np.copy(dna_hull.vertices)  # np.copy function makes an array copy of object
    dna_hull_points = np.copy(points[dna_vertices])  # np.copy function makes an array copy of object

    # Protein hull
    prot_hull = ConvexHull(prot_points)
    prot_vertices = np.copy(prot_hull.vertices)  # np.copy function makes an array copy of object
    prot_hull_points = np.copy(points[prot_vertices])  # np.copy function makes an array copy of object

    # Intersection hull


    # Find common points between DNA hull and protein hull points
    common_points = []

    for prot_hull_point in prot_hull_points:
        for dna_hull_point in dna_hull_points:
            if np.all(prot_hull_point == dna_hull_point):
                common_points.append(prot_hull_point)
                break

    common_points = np.array(common_points)

    # Create a ConvexHull object for the intersection hull of common points
    intersection_hull = ConvexHull(common_points)

    # Extract the intersection hull points
    intersection_hull_points = common_points[intersection_hull.vertices]
    # Create figure
    fig = go.Figure()
    # Add scatter points with opacity
    fig.add_trace(
        go.Scatter3d(x=points[:, 0], y=points[:, 1], z=points[:, 2], mode='markers', marker=dict(size=3, color='gray'),
                     name='All Points', opacity=0.125))
    fig.add_trace(go.Scatter3d(x=prot_points[:, 0], y=prot_points[:, 1], z=prot_points[:, 2], mode='lines',
                               line=dict(color='red'), name='Transcription Factor'))
    fig.add_trace(
        go.Scatter3d(x=dna_points[:, 0], y=dna_points[:, 1], z=dna_points[:, 2], mode='lines', line=dict(color='blue'),
                     name='Nucleic Acid'))

    # Add convex hulls with opacity
    #fig.add_trace(
        #go.Mesh3d(x=dna_hull_points[:, 0], y=dna_hull_points[:, 1], z=dna_hull_points[:, 2], color='blue', opacity=0.25,
                 # name='DNA Hull'))
    #fig.add_trace(go.Mesh3d(x=prot_hull_points[:, 0], y=prot_hull_points[:, 1], z=prot_hull_points[:, 2], color='red',
                            #opacity=0.25, name='Protein Hull'))
    #fig.add_trace
    # (x=intersection_hull_points[:, 0], y=intersection_hull_points[:, 1], z=intersection_hull_points[:, 2],
                  #color='purple', opacity=0.5, name='Interaction Site Hull'))
    fig.add_trace(
        go.Mesh3d(x=hull_points[:, 0], y=hull_points[:, 1], z=hull_points[:, 2],
                  color='purple', opacity=0.5, name='Complex Hull'))


    # Set layout
    fig.update_layout(scene=dict(aspectmode="data"))
    fig.show()