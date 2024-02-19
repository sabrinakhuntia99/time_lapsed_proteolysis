import plotly.graph_objs as go
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

def plot_hulls_for_time_points(conv_peptide_coord, hull_layers):
    # Initialize lists to store time points and corresponding hull numbers detected
    time_points = []
    hull_numbers = []

    # Iterate over each data point in the conv_peptide_coord data
    for data_point in conv_peptide_coord:
        time_point = data_point['TimePoint']
        point = (data_point['X'], data_point['Y'], data_point['Z'])

        # Find the hull that the peptide point belongs to
        found_hull = None
        for j, hull_points in enumerate(hull_layers):
            hull = ConvexHull(hull_points)
            all_on_same_side = all(
                np.dot(eq[:-1], point) + eq[-1] <= 0 for eq in hull.equations
            )
            if all_on_same_side:
                found_hull = j + 1
                break

        # Append the time point and hull number to the lists
        time_points.append(time_point)
        hull_numbers.append(found_hull)

    # Plotting the scatter graph
    plt.figure(figsize=(8, 6))
    plt.scatter(time_points, hull_numbers, marker='o', color='blue')
    plt.xlabel('Time Points')
    plt.ylabel('Detected Hull Number')
    plt.title('Detected Hull Numbers for Peptides over Time')

    # Calculate the maximum number of hull layers
    max_hull_number = len(hull_layers)

    # Set the y-axis limit in increments of 1 up to the maximum hull layer number
    plt.ylim(0, max_hull_number + 1)

    # Set y-axis ticks in increments of 1
    plt.yticks(range(1, max_hull_number + 2))

    plt.grid(True)
    plt.show()

def plot_xyz_coordinates_with_time_and_hulls(data, hull_layers):
    # Create an empty list to store traces
    traces = []

    # Create an empty dictionary to store detected points based on time point
    detected_points = {}

    # Iterate over the data list to group detected points based on time point
    for time_point_data in data:
        # Extract time point name and coordinates
        time_point = time_point_data['TimePoint']
        X = time_point_data['X']
        Y = time_point_data['Y']
        Z = time_point_data['Z']

        # Add coordinates to the detected points dictionary
        if time_point in detected_points:
            detected_points[time_point].append((X, Y, Z))
        else:
            detected_points[time_point] = [(X, Y, Z)]

    # Create traces for detected points and segments for each time point
    for time_point, points in detected_points.items():
        # Sort points based on X coordinate (assuming X represents the detected direction)
        points.sort(key=lambda p: p[0])

        # Extract X, Y, Z coordinates from sorted points
        X, Y, Z = zip(*points)

        # Create scatter trace for detected points
        trace_points = go.Scatter3d(
            x=X,
            y=Y,
            z=Z,
            mode='markers',
            marker=dict(size=8, color='white', symbol='square'),  # Adjust marker size and color as needed
            name=f'Points at Time Point: {time_point}'  # Include time point name in trace name
        )
        traces.append(trace_points)

        # Create line trace for detected segments
        trace_detected = go.Scatter3d(
            x=X,
            y=Y,
            z=Z,
            mode='lines',
            line=dict(color='white', width=4),  # Adjust line color and width as needed
            name=f'Peptide at Time Point: {time_point}'  # Include time point name in trace name
        )
        traces.append(trace_detected)

    # Define a list of colors for the hulls
    hull_colors = ['red', 'green', 'blue', 'yellow', 'orange', 'purple', 'cyan', 'magenta', 'lime', 'pink']

    # Create a trace for each hull layer with a unique color
    for i, hull_points in enumerate(hull_layers):
        # Use modulo operation to repeat colors if there are more hulls than colors
        color_index = i % len(hull_colors)
        color = hull_colors[color_index]

        trace_hull = go.Mesh3d(
            x=hull_points[:, 0],
            y=hull_points[:, 1],
            z=hull_points[:, 2],
            opacity=0.3,  # Adjust opacity as needed
            color=color,  # Assign a unique color to each hull
            name=f'Hull {i + 1}'
        )
        traces.append(trace_hull)

    # Create plot layout

    layout = go.Layout(
        title='Convex Layered Protein Structure with Detected Peptides (Cartesian Coordinates)',
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

    # Create the figure
    fig = go.Figure(data=traces, layout=layout)

    # Show the plot
    fig.show()
