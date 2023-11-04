# Mathematical Techniques for Structural Proteomics
## Applying Convex Hull Layering to Time-Lapsed Proteolysis

convex hull = smallest convex shape that contains the given set of points

Let S be a set of points (for biomolecules, S = atomic coordinates extracted from a .pdb file).
Let P be a subset of points (for proteins, P = nitrogen coordinates, AKA the backbone, extracted from a .pdb file).

### To find a 3D convex hull (Python has a package to compute it):

1. Determine the largest convex polyhedron whose vertices (most extreme observations, by definition) are points in P.
2. Determine the set of all convex combinations of points in P.
3. Determine the intersection of all convex regions containing P.

### To continue peeling:
#### Nested Series Method
4. "Delete" vertices (for calculation purposes only -- this step accounts for the extreme observations).
5. Repeat steps 1-3 (in other words, keep finding the new hull based on the previous hull) until the minimum hull size is reached.

#### Radial Project Method
4. Compute a matrix of distances between the centroid and each data point.
5. Identify gaps as outer boundaries of shells.

### To analyze the enzymatic mechanism:
6. Designate which layers the cleavage sites are within at each checkpoint.

### To visualize results:
7. Plot:
- Atomic coordinates in gray
- Backbone coordinates in black
- Convex hull layers in multiple colors (translucent)
- Cleaved peptides in red
   
