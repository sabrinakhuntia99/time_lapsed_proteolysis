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