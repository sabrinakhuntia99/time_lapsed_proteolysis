from scipy.spatial import ConvexHull
import numpy as np
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
            layers.append(hull_points)
        except:
            break

    return layers

def calculate_average_depth_of_hulls(layers):
    total_depth = 0
    num_layers = len(layers)

    for layer in layers:
        # Calculate the depth of each layer as the difference between the maximum and minimum z-coordinate
        max_z = np.max(layer[:, 2])
        min_z = np.min(layer[:, 2])
        depth = max_z - min_z
        total_depth += depth

    # Calculate the average depth
    if num_layers > 0:
        average_depth = total_depth / num_layers
    else:
        average_depth = 0

    return average_depth
