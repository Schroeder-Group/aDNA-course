import pandas as pd

# Load your data
# Replace '/path/to/your/file.csv' with the path to your CSV file
data = pd.read_csv('319publish_Lola_lepromatosis.csv')

# Open a file to write the FASTA output
with open('319publish_Lola_lepromatosis.fasta', 'w') as fasta_file:
    # Iterate through each row of the dataframe
    for index, row in data.iterrows():
        # The first column is the sample name, which we use in the header
        header = f">{row[0]}"
        # Concatenate all remaining columns to form the sequence
        sequence = ''.join(row[1:])
        # Write to the FASTA file
        fasta_file.write(f"{header}\n{sequence}\n")

print("FASTA file has been created successfully.")

