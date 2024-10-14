import pandas as pd
from Bio import SeqIO
import sys

# Function to process NetMHC txt output with multiple header rows
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
    
    # Print the final el_rank_columns for verification
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

    return selected_data

def parse_fasta_file(fasta_file_path):
    fasta_data = []
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        header = record.id
        sequence = str(record.seq)
        parts = header.split(':')
        gene_name = parts[0]
        gene_id = parts[1]
        nchange = parts[3]
        achange = parts[4].split('_')[0]
        status = 'mut' if 'mut' in parts[4] else 'wt'
        fasta_data.append({
            'Gene': gene_name,
            'GeneID': gene_id,
            'Nchange': nchange,
            'Achange': achange,
            'Epi9mer': sequence,
            'status': status
        })
    return pd.DataFrame(fasta_data)

def merge_data(netmhc_data, fasta_data):
    return pd.merge(fasta_data, netmhc_data, left_on='Epi9mer', right_on='Peptide', how='left')

def main(netmhc_file_path, fasta_file_path, output_file_path):
    netmhc_data = process_netmhc_txt(netmhc_file_path)
    fasta_data = parse_fasta_file(fasta_file_path)

    merged_data = merge_data(netmhc_data, fasta_data)

    final_data = merged_data.drop_duplicates()

    final_data.to_excel(output_file_path, index=False)
    print(f"\nFinal formatted data has been written to {output_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python format.py [netmhc_file_path] [fasta_file_path] [output_file_path]")
        sys.exit(1)

    netmhc_file_path = sys.argv[1]
    fasta_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    main(netmhc_file_path, fasta_file_path, output_file_path)
