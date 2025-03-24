import pandas as pd
import numpy as np
from Bio import SeqIO
import sys
import re
from io import StringIO

# --------------------------------------------------
# Function to process NetMHC txt output with multiple header rows
# --------------------------------------------------
def process_netmhc_txt(file_path):
    # Read the first row separately to get the HLA alleles
    hla_alleles = pd.read_csv(file_path, delimiter='\t', nrows=1, header=None).iloc[0]

    # Clean up the HLA alleles, keeping only valid ones
    hla_alleles = hla_alleles[hla_alleles.str.startswith("HLA", na=False)]
    
    # Re-number the indices of hla_alleles after cleanup
    hla_alleles = hla_alleles.reset_index(drop=True)
    
    # Read the actual data starting from the second row
    raw_data = pd.read_csv(file_path, delimiter='\t', skiprows=1)
    
    # Print out the actual column names after loading the data
    print("Actual column names in raw data:")
    print(raw_data.columns.tolist())

    # Print the cleaned hla_alleles for debugging
    print("\nHLA Alleles:")
    print(hla_alleles)

    # Initialize a list for renaming the EL_Rank columns
    el_rank_columns = []
    new_column_names = ['Peptide']  # Start with required columns

    # Iterate through each column in the raw data and match with HLA alleles
    hla_idx = 0
    for i, col in enumerate(raw_data.columns):
        if 'EL_Rank' in col:
            if hla_idx < len(hla_alleles):
                hla = hla_alleles[hla_idx]  # Get the corresponding HLA allele from the cleaned hla_alleles
                new_name = f'ELrank_{hla.replace("-", "_")}'
                new_column_names.append(new_name)
                el_rank_columns.append(col)
                # Print statements for debugging
                print(f"Iteration {i}:")
                print(f"HLA: {hla}, New Name: {new_name}")
                hla_idx += 1
    
    # Print the final EL_Rank Columns for verification
    print(f"\nFinal EL_Rank Columns: {el_rank_columns}")
    # Print the final selected columns for verification
    print(f"\nFinal selected columns: {new_column_names}")

    # Select the necessary columns
    selected_data = raw_data[['Peptide'] + el_rank_columns]
    # Rename the selected EL_Rank columns
    selected_data.columns = new_column_names

    # Print the first 5 rows for troubleshooting
    print("\nFirst 5 rows of the processed NetMHC data:")
    print(selected_data.head())

    # For our merge strategy, add a "block" column.
    # Since each original FASTA peptide produces 26 rows, we set:
    block_size = 26
    selected_data.loc[:, 'block'] = np.arange(len(selected_data)) // block_size

    return selected_data

# --------------------------------------------------
# Function to parse FASTA file with mutation information
# --------------------------------------------------
def parse_fasta_file(fasta_file_path):
    fasta_data = []
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        header = record.id
        sequence = str(record.seq)
        parts = header.split(':')
        gene_name = parts[0]
        gene_id = parts[1]
        nchange = parts[3] if len(parts) > 3 else ""
        # parts[4] may contain additional info separated by underscores.
        if len(parts) > 4:
            achange = parts[4].split('_')[0]
        else:
            achange = ""
        # For status, check if the header contains 'mut' or 'wt'; if not, default to 'mut'
        if 'mut' in header.lower():
            status = 'mut'
        elif 'wt' in header.lower():
            status = 'wt'
        else:
            status = 'mut'
        fasta_data.append({
            'Gene': gene_name,
            'GeneID': gene_id,
            'Nchange': nchange,
            'Achange': achange,
            'Epi15mer': sequence,
            'status': status
        })
    fasta_df = pd.DataFrame(fasta_data)
    # Assign block number corresponding to the FASTA record order (starting at 0)
    fasta_df['block'] = fasta_df.index
    return fasta_df

# --------------------------------------------------
# Function to merge NetMHC output and FASTA data using block
# --------------------------------------------------
def merge_data(netmhc_df, fasta_df):
    merged_df = pd.merge(fasta_df, netmhc_df, on='block', how='left', suffixes=('_fasta', '_netMHC'))
    return merged_df

# --------------------------------------------------
# Main processing function
# --------------------------------------------------
def main(netmhc_file_path, fasta_file_path, output_file_path):
    netmhc_data = process_netmhc_txt(netmhc_file_path)
    fasta_data = parse_fasta_file(fasta_file_path)

    merged_data = merge_data(netmhc_data, fasta_data)

    # Write the merged data to a CSV file (tab-separated) instead of Excel to avoid sheet size issues.
    csv_output = output_file_path.rsplit('.', 1)[0] + ".csv"
    merged_data.to_csv(csv_output, index=False, sep='\t')
    print(f"\nFinal formatted data has been written to {csv_output}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python format_merge.py [netmhc_txt_file] [fasta_file] [output_file]")
        sys.exit(1)

    netmhc_file_path = sys.argv[1]
    fasta_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    main(netmhc_file_path, fasta_file_path, output_file_path)
