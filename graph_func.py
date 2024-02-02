# Initialize a list to store line segments as tuples of start and end points
peptide_lines = []

# Extract coordinates for peptide O60361 from the peptideList
for peptide in peptideList:
    if peptide[0] == 'O60361':  # Check for the peptide name
        peptide_O60361_available = True  # Set the flag if found
        for i in peptide:
            if i != 'name':  # Skip the 'name' key
                pointsInPeptide = peptide[i]
                for point in pointsInPeptide
            peptide_lines.extend([start_point, end_point])
        print(peptide_lines)

if not peptide_available:
    print("Peptide data not found.")
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