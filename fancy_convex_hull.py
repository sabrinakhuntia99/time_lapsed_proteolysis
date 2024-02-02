from __future__ import division
from scipy.spatial import ConvexHull
from Bio.PDB import *
import numpy as np
import csv
import plotly.graph_objs as go
import ast
import re

p = PDBParser()

# AlphaFold structure for tests
structure = p.get_structure('test_1', r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\test_1.pdb")
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

# Merge DNA and protein points
combined_points = np.concatenate((dna_points, prot_points))

# Compute convex hull for combined points
hull = ConvexHull(combined_points)
vertices = np.copy(hull.vertices)
hull_points = np.copy(combined_points[vertices])

# Create a Plotly 3D scatter plot for combined DNA and protein coordinates with blue DNA portion
trace_combined = go.Scatter3d(
    x=combined_points[:, 0],
    y=combined_points[:, 1],
    z=combined_points[:, 2],
    mode='markers+lines',  # Use markers and lines for the combined coordinates
    marker=dict(size=5, color='rgba(0, 0, 255, 0.8)', opacity=0.8, symbol='circle'),  # Blue circles for DNA portion
    line=dict(color='rgba(0, 0, 255, 0.6)', width=2),  # Translucent lines for connecting DNA points
    name='Combined Coordinates'
)

# Create Plotly 3D surface plot for the single hull with purple color
trace_hull = go.Mesh3d(
    x=hull_points[:, 0],
    y=hull_points[:, 1],
    z=hull_points[:, 2],
    opacity=0.2,  # Set lower opacity for the hull
    color='rgba(128, 0, 128, 0.6)',  # Purple color for the hull
    lighting=dict(ambient=0.6, diffuse=0.8),  # Adjust lighting properties
    name='Combined Hull'
)

layout = go.Layout(
    scene=dict(
        xaxis=dict(title='X', backgroundcolor="black", gridcolor="gray", showgrid=False),
        yaxis=dict(title='Y', backgroundcolor="black", gridcolor="gray", showgrid=False),
        zaxis=dict(title='Z', backgroundcolor="black", gridcolor="gray", showgrid=False),
        bgcolor='black',
        aspectmode="cube",
    ),
    paper_bgcolor='black',
    plot_bgcolor='black',
    showlegend=True,
    legend=dict(x=0.85, y=0.95),
)

fig = go.Figure(data=[trace_combined, trace_hull], layout=layout)
fig.update_layout(
    title='Fancy 3D Plot',
    scene=dict(
        xaxis=dict(title='', tickfont=dict(size=10, color='white'), showgrid=False, showticklabels=False),
        yaxis=dict(title='', tickfont=dict(size=10, color='white'), showgrid=False, showticklabels=False),
        zaxis=dict(title='', tickfont=dict(size=10, color='white'), showgrid=False, showticklabels=False),
        bgcolor='black',
        camera=dict(eye=dict(x=1.2, y=1.2, z=0.6))
    ),
    legend=dict(x=0.1, y=0.9, font=dict(size=10, color='lightgrey')),
    margin=dict(l=0, r=0, b=0, t=40),
    plot_bgcolor='black',
    paper_bgcolor='black'
)

fig.update_layout(legend=dict(bgcolor='black', bordercolor='black', borderwidth=0))
fig.show()
