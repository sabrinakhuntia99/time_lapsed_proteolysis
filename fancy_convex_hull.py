import numpy as np
import plotly.graph_objs as go
from scipy.spatial import ConvexHull
from Bio.PDB import *

p = PDBParser()

# AlphaFold structure for tests
structure = p.get_structure('O60361',
                            r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\AF-O60361-F1-model_v4.pdb")

prot_points = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "N":
                    prot_point = atom.get_coord()
                    prot_points.append(prot_point)
prot_points = np.array(prot_points)

# Compute convex hull for protein points
hull = ConvexHull(prot_points)
vertices = np.copy(hull.vertices)
hull_points = np.copy(prot_points[vertices])

# Duplicate first point to ensure a closed surface
hull_points = np.vstack([hull_points, hull_points[0]])

# Create triangular faces for the hull surface
hull_faces = []
for i in range(len(hull.vertices) - 1):
    hull_faces.append([i, i + 1, len(hull.vertices)])
hull_faces.append([len(hull.vertices) - 1, 0, len(hull.vertices)])

# Create a Plotly 3D scatter plot for protein coordinates
trace_prot = go.Scatter3d(
    x=prot_points[:, 0],
    y=prot_points[:, 1],
    z=prot_points[:, 2],
    mode='markers+lines',  # Use markers and lines for the prot coordinates
    marker=dict(size=5, color='rgba(0, 0, 255, 0.8)', opacity=0.8, symbol='circle'),  # Blue circles
    line=dict(color='rgba(0, 0, 255, 0.6)', width=2),  # Translucent lines
    name='Protein Coordinates'
)

# Create a closed surface trace for the convex hull
trace_hull = go.Mesh3d(
    x=hull_points[:, 0],
    y=hull_points[:, 1],
    z=hull_points[:, 2],
    i=hull_faces,
    j=hull_faces,
    k=hull_faces,
    opacity=0.3,  # Adjust opacity as needed
    color='blue',  # Assign a unique color to the hull
    name='Convex Hull'
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

fig = go.Figure(data=[trace_prot, trace_hull], layout=layout)
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
