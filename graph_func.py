import plotly.graph_objs as go

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

