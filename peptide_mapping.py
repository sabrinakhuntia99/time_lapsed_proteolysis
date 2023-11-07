import pandas as pd
import ahocorasick as acs
from bokeh.plotting import figure, show, output_file
import gzip
import csv

# extract peptide distance values
peptides = []
with open("peptides.tsv") as file:
    reader = csv.reader(file, delimiter="\t")
    for row in reader:
        peptides.append(row)

peptides_df = peptides[:0]

# Aho-Corasick automaton for peptides
peptides = acs.Automaton()
for peptide in peptides_df['Peptide']:
    peptides.add_word(peptide, peptide)
peptides.make_automaton()

# dictionary to store peptide-protein mappings
matched_peptides = {}

# Read in proteins from UniProt_Human.fasta.gz
with gzip.open(r"UniProt_Human.fasta.gz", "rt") as fin:
    protein_sequences = fin.read()

    # Iterate over matches
    for end_index, peptide in peptides.iter(protein_sequences):
        start_index = end_index - len(peptide) + 1

        # Extract protein entry from protein sequences
        protein_entry = protein_sequences.rfind('>', 0, start_index)
        protein_entry = protein_sequences[protein_entry:end_index].strip()

        # Extract protein ID and sequence
        lines = protein_entry.split('\n')
        protein_id = lines[0].split()[0]
        protein_sequence = ''.join(lines[1:])

        if peptide in protein_sequence:
            if peptide not in matched_peptides:
                matched_peptides[peptide] = []
            matched_peptides[peptide].append(protein_id)

# Print peptide-protein mappings
for peptide, proteins in matched_peptides.items():
    print(f"Peptide: {peptide}, Proteins: {', '.join(proteins)}")


# Read the peptide column from psm.tsv.gz
with gzip.open(r"psm.tsv.gz", "rt") as fin:
    peptides_df = pd.read_csv(fin, sep='\t', usecols=['Peptide'])

# Create Aho-Corasick automaton for peptides
peptides = acs.Automaton()
for peptide in peptides_df['Peptide']:
    peptides.add_word(peptide, peptide)
peptides.make_automaton()

# Create dictionary to store peptide-protein mappings
matched_peptides = {}

# Read in proteins from UniProt_Human.fasta.gz
with gzip.open(r"UniProt_Human.fasta.gz", "rt") as fin:
    protein_sequences = fin.read()

    # Iterate over the matches found by the automaton
    for end_index, peptide in peptides.iter(protein_sequences):
        start_index = end_index - len(peptide) + 1

        # Extract protein entry from the protein sequences
        protein_entry = protein_sequences.rfind('>', 0, start_index)
        protein_entry = protein_sequences[protein_entry:end_index].strip()

        # Extract protein ID and sequence
        lines = protein_entry.split('\n')
        protein_id = lines[0].split()[0]
        protein_sequence = ''.join(lines[1:])

        if peptide in protein_sequence:
            if peptide not in matched_peptides:
                matched_peptides[peptide] = []
            matched_peptides[peptide].append(protein_id)

# Print peptide-protein mappings
for peptide, proteins in matched_peptides.items():
   print(f"Peptide: {peptide}, Proteins: {', '.join(proteins)}")

# frequency of each identified protein
protein_counts = {}
for proteins in matched_peptides.values():
    for protein in proteins:
        if protein not in protein_counts:
            protein_counts[protein] = 0
        protein_counts[protein] += 1

# dataframe for the protein counts
protein_df = pd.DataFrame(list(protein_counts.items()), columns=['Protein', 'Frequency'])

# Sort the dataframe by frequency in descending order
protein_df = protein_df.sort_values(by='Frequency', ascending=False)

# Bokeh figure
p = figure(x_range=protein_df['Protein'], height=600, width=1000, title="Protein Frequency Coverage Map",
           toolbar_location=None, tools="")

# appearance of the bar graph
p.vbar(x='Protein', top='Frequency', width=0.9, source=protein_df, line_color='white')

# rotate x-axis labels
p.xaxis.major_label_orientation = 1.2

# Set axis labels and plot title
p.xaxis.axis_label = "Protein ID"
p.yaxis.axis_label = "Frequency"
p.title.align = "center"

# Show the plot
show(p)
