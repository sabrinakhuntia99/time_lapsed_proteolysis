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

    # Iterate over the data list
    for time_point_data in data:
        # Extract time point name and coordinates
        time_point = time_point_data['TimePoint']
        X = time_point_data['X']
        Y = time_point_data['Y']
        Z = time_point_data['Z']

        # Create a scatter trace for the coordinates
        trace = go.Scatter3d(
            x=[X],
            y=[Y],
            z=[Z],
            mode='markers',
            marker=dict(size=5, color='blue'),  # Adjust marker size and color as needed
            name=f'Time Point: {time_point}'  # Include the time point name in the trace name
        )
        traces.append(trace)

    # Create a trace for each hull layer
    for i, hull_points in enumerate(hull_layers):
        trace_hull = go.Mesh3d(
            x=hull_points[:, 0],
            y=hull_points[:, 1],
            z=hull_points[:, 2],
            opacity=0.3,  # Adjust opacity as needed
            color='red',  # Adjust hull color as needed
            name=f'Hull {i + 1}'
        )
        traces.append(trace_hull)

    # Create layout for the plot
    layout = go.Layout(
        title='XYZ Coordinates with Time and Hulls',
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Z')
        )
    )

    # Create the figure
    fig = go.Figure(data=traces, layout=layout)

    # Show the plot
    fig.show()

