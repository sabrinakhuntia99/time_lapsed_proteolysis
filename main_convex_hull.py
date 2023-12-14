
from __future__ import division
from scipy.spatial import ConvexHull, Delaunay
from Bio.PDB import *
import numpy as np
import csv
import plotly.graph_objects as go
import matplotlib.pyplot as plt


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


def convex_hull_peeling(points):
    hull = ConvexHull(points)
    vertices = np.copy(hull.vertices)
    hull_points = np.copy(points[vertices])
    layers = [hull_points]

    min_points_percentage = 0.1  # Adjust this percentage as needed
    min_points = int(len(points) * min_points_percentage)

    while len(hull_points) > min_points:
        try:
            new_hull = ConvexHull(hull_points)
            vertices = np.copy(new_hull.vertices)
            hull_points = np.copy(points[vertices])
            print(len(hull_points))
            layers.append(hull_points)
        except:
            break

    return layers

def calculate_max_distance_to_hull(centroid, hull_points):
    return max(np.linalg.norm(p - centroid) for p in hull_points)

def check_vectors_within_max_distance(vectors, max_distance):
    within_max_distance = [np.linalg.norm(vector) <= max_distance for vector in vectors]
    return within_max_distance



if __name__ == "__main__":

    # Perform convex hull peeling
    layers = convex_hull_peeling(points)
    print("Number of Hulls:", len(layers))

    # Initialize lists to store coordinates
    x_coordinates = []
    y_coordinates = []
    z_coordinates = []

    # Initialize lists to store line segment coordinates
    peptide_lines_x = []
    peptide_lines_y = []
    peptide_lines_z = []

    # Extract coordinates for peptide O60361 from the peptideList
    peptide_coordinates = []  # List to store extracted coordinates

    # Initialize a list to store line segments as tuples of start and end points
    peptide_lines = []

    # Extract coordinates for peptide O60361 from the peptideList
    for peptide in peptideList:
        if peptide['name'] == 'O60361':  # Check for the peptide name
            peptide_O60361_available = True  # Set the flag if found
            for i in peptide:
                if i != 'name':  # Skip the 'name' key
                    pointsInPeptide = peptide[i]
                    for point in pointsInPeptide:
                        currentPoint = pointsInPeptide[point]
                        # Extract coordinates and append to the list
                        peptide_coordinates.append((currentPoint['x'], currentPoint['y'], currentPoint['z']))

            # Create line segments from the extracted coordinates
            for i in range(0, len(peptide_coordinates) - 1, 2):
                start_point = peptide_coordinates[i]
                end_point = peptide_coordinates[i + 1]

                # Append each start and end point as a separate tuple
                peptide_lines.extend([start_point, end_point])
            print(peptide_lines)

    if not peptide_O60361_available:
        print("Peptide O60361 data not found.")
    else:
        # Create traces for line segments of peptide O60361
        trace_peptide_lines = go.Scatter3d(
            x=[point[0] for point in peptide_lines],  # X-coordinates for start and end points
            y=[point[1] for point in peptide_lines],  # Y-coordinates for start and end points
            z=[point[2] for point in peptide_lines],  # Z-coordinates for start and end points
            mode='lines',
            line=dict(color='red', width=9),
            name='Detected Peptides'
        )
    # Create a trace for peptide O60361 coordinates
    trace_peptide_O60361_start = go.Scatter3d(
        x=[point[0] for point in peptide_coordinates[::2]],  # X-coordinates for start points
        y=[point[1] for point in peptide_coordinates[::2]],  # Y-coordinates for start points
        z=[point[2] for point in peptide_coordinates[::2]],  # Z-coordinates for start points
        mode='markers',
        marker=dict(size=10, color='rosybrown', opacity=1, symbol='square', line=dict(color='rosybrown', width=2)),
        name='Start Points'
    )

    trace_peptide_O60361_end = go.Scatter3d(
        x=[point[0] for point in peptide_coordinates[1::2]],  # X-coordinates for end points
        y=[point[1] for point in peptide_coordinates[1::2]],  # Y-coordinates for end points
        z=[point[2] for point in peptide_coordinates[1::2]],  # Z-coordinates for end points
        mode='markers',
        marker=dict(size=10, color='crimson', opacity=1, symbol='square', line=dict(color='crimson', width=2)),
        name='End Points'
    )

    # Printing the coordinates used in the traces
    print("Start Points Coordinates:")
    for x, y, z in zip(trace_peptide_O60361_start['x'], trace_peptide_O60361_start['y'],
                       trace_peptide_O60361_start['z']):
        print(f"({x}, {y}, {z})")

    print("\nEnd Points Coordinates:")
    for x, y, z in zip(trace_peptide_O60361_end['x'], trace_peptide_O60361_end['y'], trace_peptide_O60361_end['z']):
        print(f"({x}, {y}, {z})")


    trace_atomic = go.Scatter3d(
        x=points0[:, 0],
        y=points0[:, 1],
        z=points0[:, 2],
        mode='markers',
        marker=dict(size=3, color='grey', opacity=0.25),
        name='Atomic Coordinates'
    )

    trace_amino_acid = go.Scatter3d(
        x=points[:, 0],
        y=points[:, 1],
        z=points[:, 2],
        mode='markers',
        marker=dict(size=5, color='white', opacity=1, symbol='square', line=dict(color='white', width=2)),
        name='Amino Acid Backbone'
    )

    # Create traces for hulls
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
            colorscale=[[0, 'lightblue'], [1, 'deepskyblue']],
        )
        trace_hulls.append(trace_hull)

    # Calculate centroid of the points
    centroid = np.mean(points, axis=0)

    # Create a trace for the centroid point
    trace_centroid = go.Scatter3d(
        x=[centroid[0]],  # X-coordinate of the centroid
        y=[centroid[1]],  # Y-coordinate of the centroid
        z=[centroid[2]],  # Z-coordinate of the centroid
        mode='markers',
        marker=dict(size=10, color='purple', opacity=1, symbol='square', line=dict(color='purple', width=2)),
        name='Centroid'
    )



    # Combine relevant traces into a single list
    combined_data = [trace_atomic, trace_amino_acid] + trace_hulls + [trace_peptide_O60361_start, trace_peptide_O60361_end, trace_peptide_lines]
    combined_data_with_centroid = combined_data + [trace_centroid]

    # Create the layout
    layout = go.Layout(
        title=("UNIPROT_ID O60361"),
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

    # Create the figure with combined data
    fig = go.Figure(data=combined_data_with_centroid, layout=layout)

    # Show the interactive 3D plot
    fig.show()

    # List to store which hull each peptide point is in
    peptide_point_hull_assignments = []

    # Assuming 'peptide_coordinates' contains the points of interest within the peptide
    for point in peptide_coordinates:
        found_hull = None
        for j, hull_points in enumerate(layers):
            hull = ConvexHull(hull_points)
            all_on_same_side = all(
                np.dot(eq[:-1], point) + eq[-1] <= 0 for eq in hull.equations
            )
            if all_on_same_side:
                found_hull = j + 1
                break

        if found_hull is not None:
            peptide_point_hull_assignments.append(f"Point {point} is in Hull {found_hull}")
        else:
            peptide_point_hull_assignments.append(f"Point {point} is not in any hull")

    # Print which hulls the peptide points are detected in
    for assignment in peptide_point_hull_assignments:
        print(assignment)

    # Check which hull the centroid is in
    centroid_hull_assignment = None
    for i, hull_points in enumerate(layers):
        hull = ConvexHull(hull_points)
        is_inside_hull = all(
            np.dot(eq[:-1], centroid) + eq[-1] <= 0 for eq in hull.equations
        )
        if is_inside_hull:
            centroid_hull_assignment = f"Centroid is in Hull {i + 1}"
            break

    if centroid_hull_assignment is not None:
        print(centroid_hull_assignment)
    else:
        print("Centroid is not inside any hull")

