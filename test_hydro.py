from __future__ import division
from scipy.spatial import ConvexHull, Delaunay
from Bio.PDB import *
import numpy as np
import csv
import plotly.graph_objs as go
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

# Amino Acid Hydrophobicity Scale (Kyte-Doolittle scale)
hydrophobicity_values = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    'O': 0.0, # not in this sequence
    'U': 0.0 # not in this sequence
}


p=PDBParser()

#  AlphaFold structure for test
structure=p.get_structure('AF-O60361-F1-model_v4', r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\AF-O60361-F1-model_v4.pdb")

# Extract amino acid sequence from the first model and chain
amino_acid_sequence = ""
for model in structure:
    for chain in model:
        for residue in chain:
            # Check if the residue is an amino acid (exclude water molecules, ligands, etc.)
            if residue.get_id()[0] == " " and residue.get_resname() not in ["HOH", "H2O", "WAT"]:
                amino_acid_sequence += residue.get_resname()

print("Amino Acid Sequence:", amino_acid_sequence)
# Create an array of hydrophobicity values based on the amino acid sequence
hydrophobicity_array = [hydrophobicity_values[aa] for aa in amino_acid_sequence]

print(hydrophobicity_array)

# backbone atoms
points = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "N":
                    point = atom.get_coord()
                    points.append(point)
points = np.array(points)

# all atoms
points0 = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                point0 = atom.get_coord()
                points0.append(point0)
points0 = np.array(points0)

# extract peptide distance values
peptides = []
with open("peptides.tsv") as file:
    reader = csv.reader(file, delimiter="\t")
    for row in reader:
        peptides.append(row)

peptide = peptides[2]
vectors = peptide[1:13]
# print(f"Time-Lapsed Peptide-to-Centroid Distances: {vectors}")

input_line = []

with open("pep_cleave_coordinates_10292023.csv") as file:
    reader = csv.reader(file, delimiter=",")
    for row in reader:
        input_line.append(row)

peptideList = []
columns = input_line[0]
for line in input_line[1::]:
    indexCounter = 0
    peptideDict = {}
    for index in line:
        pointsDict = {}
        if (indexCounter == 0):
            peptideDict["name"] = index
        else:
            if index != '':
                if (")," in index):
                    pts = index.split("),")
                    startingPoints = pts[0]
                    startingPoints = startingPoints[8:len(startingPoints) - 1]
                    startingPoints = startingPoints.split(",")
                    xyzDict = {}
                    xyzDict["x"] = float(startingPoints[0])
                    xyzDict["y"] = float(startingPoints[1])
                    xyzDict["z"] = float(startingPoints[2])
                    pointsDict["Points1"] = xyzDict
                    pointsCounter = 2
                    for point in pts[1:len(pts) - 1]:
                        middlePoints = point[8:len(point) - 1]
                        middlePoints = middlePoints.split(",")
                        xyzDict = {}
                        xyzDict["x"] = float(middlePoints[0])
                        xyzDict["y"] = float(middlePoints[1])
                        xyzDict["z"] = float(middlePoints[2])
                        indexName = "Points" + str(pointsCounter)
                        pointsCounter += 1
                        pointsDict[indexName] = xyzDict

                    endingPoints = pts[len(pts) - 1]
                    endingPoints = endingPoints[8:len(endingPoints) - 3]
                    endingPoints = endingPoints.split(",")
                    xyzDict = {}
                    xyzDict["x"] = float(endingPoints[0])
                    xyzDict["y"] = float(endingPoints[1])
                    xyzDict["z"] = float(endingPoints[2])
                    indexName = "Points" + str(pointsCounter)
                    pointsDict[indexName] = xyzDict
                else:
                    startingPoints = index[8:len(index) - 3]
                    startingPoints = startingPoints.split(",")
                    xyzDict = {}
                    xyzDict["x"] = float(startingPoints[0])
                    xyzDict["y"] = float(startingPoints[1])
                    xyzDict["z"] = float(startingPoints[2])
                    pointsDict["Points1"] = xyzDict

                peptideDict[columns[indexCounter]] = pointsDict
        indexCounter += 1
    peptideList.append(peptideDict)

xList = []
yList = []
zList = []

# for peptide in peptideList:
#     print(peptide['name'])
#     for i in peptide:
#         if i != 'name':  # Skip the 'name' key
#             pointsInPeptide = peptide[i]
#             for point in pointsInPeptide:
#                 currentPoint = pointsInPeptide[point]
#                 print('Coordinates:', currentPoint['x'], ",", currentPoint['y'], ",", currentPoint['z'])
#
# i = 0
# while i < len(xList):
#     print('Coordinates:', xList[i], ",", yList[i], ",", zList[i])
#     i += 1


peptide_O60361_available = False  # Flag to track availability

for peptide in peptideList:
    if peptide['name'] == 'O60361':  # Check for the peptide name
        peptide_O60361_available = True  # Set the flag if found
        for i in peptide:
            if i != 'name':  # Skip the 'name' key
                pointsInPeptide = peptide[i]
                for point in pointsInPeptide:
                    currentPoint = pointsInPeptide[point]
                    print('Coordinates:', currentPoint['x'], ",", currentPoint['y'], ",", currentPoint['z'])

if not peptide_O60361_available:
    print("Peptide O60361 data not found.")


#print(peptideList)

def convex_hull_peeling(points):
    #print(len(points))
    hull = ConvexHull(points)
    vertices = np.copy(hull.vertices) # np.copy function makes an array copy of object
    hull_points = np.copy(points[vertices]) # np.copy function makes an array copy of object
    #print(len(hull_points)) # value should be <= len(points)
    layers = [hull_points] # defines list
    # change to
    # while distance b/w any 2 points w/in the hull >= smallest length of amino acid in structure
    while len(hull_points) > 11:
        try:
            new_hull = ConvexHull(hull_points) # calculate convex hull
            vertices =  np.copy(new_hull.vertices) # makes an array copy of new hull vertices
            hull_points = np.copy(points[vertices]) # makes an array copy of new hull coordinates
            print(len(hull_points))
            layers.append(hull_points) # add to previous layer
        except:
            break

    return layers

def calculate_max_distance_to_hull(centroid, hull_points):
    return max(np.linalg.norm(p - centroid) for p in hull_points)

def check_vectors_within_max_distance(vectors, max_distance):
    within_max_distance = [np.linalg.norm(vector) <= max_distance for vector in vectors]
    return within_max_distance

# Function for Delaunay tessellation peeling

def delaunay_tessellation_peeling(points):
    tri = Delaunay(points)
    simplices = np.copy(tri.simplices)
    shells = [points[simplices]]
    while len(simplices) > 4:
        try:
            tri = Delaunay(points[simplices.flatten()])
            simplices = np.copy(tri.simplices)
            layers.append(points[simplices])
        except:
            break
    return shells





if __name__ == "__main__":

    # Apply Delaunay triangulation
    tri = Delaunay(points)

    # Get the vertices of the tetrahedra
    tetrahedra = points[tri.simplices]

    # Create a 3D plot using Matplotlib
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the Delaunay triangulation in Matplotlib
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='red', marker='o', s=50)
    ax.add_collection3d(Poly3DCollection(tetrahedra, facecolors='cyan', linewidths=1, edgecolors='r', alpha=0.25))


    ax.set_title('Delaunay Triangulation')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    ax.grid()

    plt.show()

    num_points = len(points)

    # Perform convex hull peeling)
    layers = convex_hull_peeling(points)
    # Define a uniform thickness value for all hulls
    # uniform_thickness = 0.5

    # Perform Delaunay tessellation peeling
    delaunay_layers = delaunay_tessellation_peeling(points)

    # Check if each vector is within the maximum distance for each layer
    for vector in vectors:
        for layer in layers:
            # Calculate the centroid of the current hull layer
            centroid = np.mean(layer, axis=0)

            # Calculate the sets of boundary points for the current hull
            hull = ConvexHull(layer)
            boundary_points = layer[hull.vertices]

            # Calculate the maximum distance from the centroid to any point on the hull
            max_distance = calculate_max_distance_to_hull(centroid, boundary_points)  # Define max_distance here
            # adjusted_max_distance = max_distance + uniform_thickness

            # Check if each vector is within the maximum distance
            within_max_distance = check_vectors_within_max_distance(vectors, max_distance)
            # adj_within_max_distance = check_vectors_within_max_distance(vectors, adjusted_max_distance)

            # Print the results for the current hull layer
            #print(f"Hull Boundary (from Centroid): {max_distance}")
            #print(f"Vector Distances (from Centroid): {vectors}")
            #print(f"Peptide Detection within Hull by Timepoint: {within_max_distance}")

    # Extract XYZ coordinates from peptideList
    peptide_x = []
    peptide_y = []
    peptide_z = []

    for peptide_dict in peptideList:
        for key, value in peptide_dict.items():
            if key.startswith('pts'):
                peptide_x.append(value['x'])
                peptide_y.append(value['y'])
                peptide_z.append(value['z'])

    # Convert the peptide coordinates to NumPy arrays
    peptide_x = np.array(peptide_x)
    peptide_y = np.array(peptide_y)
    peptide_z = np.array(peptide_z)

    # Create a Plotly  plot for peptide coordinates
    trace_peptide = go.Scatter3d(
        x=peptide_x,
        y=peptide_y,
        z=peptide_z,
        mode='markers',
        marker=dict(size=3, color='green', opacity=0.8),
        name='Peptide Coordinates'
    )

    # Extract coordinates for peptide O60361
    peptide_O60361_x = []
    peptide_O60361_y = []
    peptide_O60361_z = []

    # Extract coordinates for peptide O60361 from the peptideList
    for peptide_dict in peptideList:
        if peptide_dict['name'] == 'O60361':  # Check for the peptide name
            for key, value in peptide_dict.items():
                if key.startswith('pts'):
                    peptide_O60361_x.append(value['x'])
                    peptide_O60361_y.append(value['y'])
                    peptide_O60361_z.append(value['z'])

    # Convert the peptide O60361 coordinates to NumPy arrays
    peptide_O60361_x = np.array(peptide_O60361_x)
    peptide_O60361_y = np.array(peptide_O60361_y)
    peptide_O60361_z = np.array(peptide_O60361_z)
    print('O60361 Coordinates:', peptide_O60361_x, ",", peptide_O60361_y, ",",  peptide_O60361_z)

    # Plot peptide O60361 coordinates onto the existing convex hull plot
    trace_peptide_O60361 = go.Scatter3d(
        x=peptide_O60361_x,
        y=peptide_O60361_y,
        z=peptide_O60361_z,
        mode='markers',
        marker=dict(size=3, color='orange', opacity=0.8),
        name='Peptide O60361 Coordinates'
    )

    # Create a Plotly 3D scatter plot for atomic coordinates
    trace_atomic = go.Scatter3d(
        x=points0[:, 0],
        y=points0[:, 1],
        z=points0[:, 2],
        mode='markers',
        marker=dict(size=3, color='grey', opacity=0.25),
        name='Atomic Coordinates'
    )

    # Create a Plotly 3D line plot for amino acid coordinates
    trace_amino_acid = go.Scatter3d(
        x=points[:, 0],
        y=points[:, 1],
        z=points[:, 2],
        mode='markers',
        marker=dict(color='black', opacity=1),
        name='Amino Acid Backbone'
    )
    # Create Plotly 3D surface plots for all peeled convex hulls with different opacities and blue gradient
    trace_hulls = []
    opacities = np.linspace(0.1, 0.8, len(layers))
    for i, hull_points in enumerate(layers):
        opacity = opacities[i]
        trace_hull = go.Mesh3d(
            x=hull_points[:, 0],
            y=hull_points[:, 1],
            z=hull_points[:, 2],
            opacity=opacity,
            name=f'Hull {i + 1}',
            colorscale=[[0, 'navy'], [1, 'deepskyblue']],  # Blue gradient
        )
        trace_hulls.append(trace_hull)

    # Create a layout with scene settings for 3D interaction, setting a black background
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

    # Create the figure and add traces with square-like markers for different components
    fig = go.Figure(
        data=[trace_atomic, trace_amino_acid, trace_peptide_O60361] + trace_hulls,
        layout=layout
    )

    # Change markers to square-like for atomic coordinates
    fig.update_traces(marker=dict(symbol='square', size=3, color='grey', opacity=0.25), selector=dict(type='scatter3d'))
    # Change markers to square-like, white, and glowy for amino acid backbone
    fig.update_traces(marker=dict(symbol='square', size=5, color='white', opacity=1, line=dict(color='white', width=1)),
                      selector=dict(name='Amino Acid Backbone'))


    # Show the interactive 3D plot
    fig.show()

    # Plot peptide O60361 coordinates onto the existing convex hull plot
    trace_peptide_O60361 = go.Scatter3d(
        x=peptide_O60361_x,
        y=peptide_O60361_y,
        z=peptide_O60361_z,
        mode='markers',
        marker=dict(size=3, color='orange', opacity=0.8),
        name='Peptide O60361 Coordinates'
    )

    # Assuming you have the existing convex hull plot in the variable trace_hull
    combined_data = [trace_hull, trace_peptide_O60361]

    # Layout settings for 3D interaction
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

    # Create the figure with combined data (convex hull and peptide coordinates)
    fig = go.Figure(data=combined_data, layout=layout)

    # Show the interactive 3D plot
    fig.show()

    # Create Plotly 3D surface plots for all Delaunay layers with different opacities and blue gradient
    trace_delaunay_layers = []
    opacities = np.linspace(0.1, 0.8, len(delaunay_layers))
    for i, layer_points in enumerate(delaunay_layers):
        opacity = opacities[i]
        trace_layer = go.Mesh3d(
            x=layer_points[:, 0],
            y=layer_points[:, 1],
            z=layer_points[:, 2],
            opacity=opacity,
            name=f'Layer {i + 1}',
            colorscale=[[0, 'navy'], [1, 'deepskyblue']],  # Blue gradient
        )
        trace_delaunay_layers.append(trace_layer)

    # Create a layout with scene settings for 3D interaction, setting a black background
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

    # Create the figure and add traces with square-like markers for different components
    fig = go.Figure(
        data=[trace_atomic, trace_amino_acid, trace_peptide] + trace_delaunay_layers,  # Replace trace_hulls
        layout=layout
    )

    # Change markers to square-like for atomic coordinates
    fig.update_traces(marker=dict(symbol='square', size=3, color='grey', opacity=0.25), selector=dict(type='scatter3d'))
    # Change markers to square-like, white, and glowy for amino acid backbone
    fig.update_traces(marker=dict(symbol='square', size=5, color='white', opacity=1, line=dict(color='white', width=1)),
                      selector=dict(name='Amino Acid Backbone'))

    # Show the interactive 3D plot
    fig.show()

    # Create Plotly 3D surface plot for the entire structure (single hull)
    avg_hydrophobicity = np.mean(hydrophobicity_array)
    trace_hull_hydrophobicity = go.Mesh3d(
        x=points[:, 0],
        y=points[:, 1],
        z=points[:, 2],
        opacity=0.5,
        colorscale='Viridis',
        intensity=hydrophobicity_array,
        cmin=min(hydrophobicity_array),
        cmax=max(hydrophobicity_array),
        name='Protein Structure (Hydrophobicity Mapping)'
    )

    # Create a layout with scene settings for 3D interaction
    layout = go.Layout(
        scene=dict(
            aspectmode="data",
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Z')
        ),
        legend=dict(x=0.85, y=0.95),
        title='Protein Structure with Hydrophobicity Mapping'
    )

    # Create the figure and add traces (including the new trace_hull_hydrophobicity)
    fig = go.Figure(data=[trace_hull_hydrophobicity] + trace_hulls, layout=layout)

    # Show the interactive 3D plot
    fig.show()

# Calculate average hydrophobicity for the entire structure
avg_hydrophobicity = np.mean(hydrophobicity_array)

# Create a Plotly 3D surface plot for a single Delaunay tessellation layer
trace_delaunay_hydrophobicity = go.Mesh3d(
    x=delaunay_layers[0][:, 0],  # Assuming the first layer from Delaunay tessellation
    y=delaunay_layers[0][:, 1],
    z=delaunay_layers[0][:, 2],
    opacity=0.5,
    colorscale='Viridis',
    intensity=hydrophobicity_array,  # Use hydrophobicity values for intensity
    cmin=min(hydrophobicity_array),
    cmax=max(hydrophobicity_array),
    name='Protein Structure (Hydrophobicity Mapping)'
)

# Create a layout with scene settings for 3D interaction
layout = go.Layout(
    scene=dict(
        aspectmode="data",
        xaxis=dict(title='X'),
        yaxis=dict(title='Y'),
        zaxis=dict(title='Z')
    ),
    legend=dict(x=0.85, y=0.95),
    title='Protein Structure with Hydrophobicity Mapping'
)

# Create the figure and add traces (including the new trace_delaunay_hydrophobicity)
fig = go.Figure(data=[trace_delaunay_hydrophobicity], layout=layout)

# Show the interactive 3D plot
fig.show()

