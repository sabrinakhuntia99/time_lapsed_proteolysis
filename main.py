import requests
import pandas as pd

# Read data from TSV file
df = pd.read_csv('disprot.tsv', sep='\t')

# Print column names
print(df.columns)

# Extract values from the 'acc' column
acc_values = df['acc'].tolist()

# Printing the extracted values
print(acc_values)

for acc in acc_values:
    url = f"https://www.disprot.org/api/{acc}"
    response = requests.get(url)
    data = response.json()
    disorder_content = data.get("disorder_content")
    print(f"Disorder Content for {acc}: {disorder_content}")