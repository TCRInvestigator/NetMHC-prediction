import pandas as pd
import subprocess
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
import os
import sys
import subprocess

# Function to clean the DNA sequence by removing non-sequence lines
def clean_dna_sequence(data):
    lines = data.strip().split('\n')
    # Remove any header or annotation lines (those starting with '>')
    sequence_lines = [line for line in lines if not line.startswith('>')]
    # Join the remaining lines to form the clean sequence
    return ''.join(sequence_lines).replace(' ', '').replace('\n', '')

# Function to fetch DNA sequence using curl
def fetch_dna_sequence(url):
    try:
        result = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True)
        cleaned_sequence = clean_dna_sequence(result.stdout)
        return cleaned_sequence
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error fetching DNA sequence: {e}")

# Function to fetch protein sequence directly from NCBI
def fetch_protein_sequence(transcript_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={transcript_id}&rettype=gb&retmode=text"
    nucleotide_record = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    protein_sequence = "".join(nucleotide_record.split('/translation="')[1].split('\"')[0].split())
    return protein_sequence

# Function to apply mutation and generate epitopes
def generate_epitopes(dna_sequence, mutation_info, start_codon_index):
    # Parse mutation information
    gene, transcript, exon, dna_change, protein_change = mutation_info.split(':')
    
    # Extract the position number from dna_change, e.g., extract 6581 from c.T6581C
    mutation_position_dna = int(''.join(filter(str.isdigit, dna_change))) - 1  # Convert to 0-based index
    
    # Extract the position number from protein_change, e.g., extract 2194 from p.I2194T
    mutation_position_aa = int(''.join(filter(str.isdigit, protein_change))) - 1  # Convert to 0-based index
    
    new_aa = protein_change[-1]

    # Fetch the original protein sequence from NCBI
    original_protein_sequence = fetch_protein_sequence(transcript)

    # Apply the mutation directly to the protein sequence
    mutated_protein_sequence = (original_protein_sequence[:mutation_position_aa] + new_aa + original_protein_sequence[mutation_position_aa + 1:])

    protein_len = len(mutated_protein_sequence)

    wildtype_epitopes = []
    mutated_epitopes = []
    epitope_length = 9  # Length of the epitope

    for i in range(epitope_length):
        # Calculate the start and end for the epitope to include the mutated amino acid
        start = mutation_position_aa - i
        end = start + epitope_length

        # Ensure the epitope window is within bounds
        if start < 0:
            end -= start
            start = 0
        if end > protein_len:
            start -= end - protein_len
            end = protein_len

        wildtype_epitope_sequence = original_protein_sequence[start:end]
        mutated_epitope_sequence = mutated_protein_sequence[start:end]

        if len(wildtype_epitope_sequence) == epitope_length:  # Ensure complete epitope sequence
            wildtype_epitopes.append(wildtype_epitope_sequence)
            mutated_epitopes.append(mutated_epitope_sequence)

    return wildtype_epitopes, mutated_epitopes

# Function to write epitopes to a FASTA file
def write_epitopes_to_fasta(all_epitopes, output_file):
    with open(output_file, 'w') as f:
        for epitope in all_epitopes:
            f.write(epitope + '\n')
    print(f"Epitopes have been written to {output_file}")

# Main function to process the CSV file and generate epitopes
def process_csv_file(csv_file_path):
    # Load the CSV file containing the mutations and URLs
    df = pd.read_csv(csv_file_path)

    # List to store all epitopes for output
    all_epitopes = []

    # Process each mutation and generate epitopes
    for index, row in df.iterrows():
        aa_changes = row['AAChange'].split(',')
        for mutation_info in aa_changes:
            gene, transcript, exon, dna_change, protein_change = mutation_info.split(':')
            
            # Construct the URL for the DNA sequence
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={transcript}&rettype=fasta&retmode=text'
            
            # Fetch and clean the DNA sequence
            dna_sequence = fetch_dna_sequence(url)
            
            # Find the first ATG (start codon)
            start_codon_index = dna_sequence.find("ATG")
            if start_codon_index == -1:
                raise ValueError("No start codon (ATG) found in the DNA sequence")

            wildtype_epitopes, mutated_epitopes = generate_epitopes(dna_sequence, mutation_info, start_codon_index)
            for i, (wildtype_epitope, mutated_epitope) in enumerate(zip(wildtype_epitopes, mutated_epitopes)):
                wt_epitope_id = f"{mutation_info}_wt_epitope_{i + 1}"
                mut_epitope_id = f"{mutation_info}_mut_epitope_{i + 1}"
                all_epitopes.append(f">{wt_epitope_id}\n{wildtype_epitope}")
                all_epitopes.append(f">{mut_epitope_id}\n{mutated_epitope}")

    # Get the base name of the CSV file without the extension
    base_name = os.path.splitext(os.path.basename(csv_file_path))[0]
    output_file = f"{base_name}_epitopes_output.fasta"

    # Write all epitopes to a FASTA file
    write_epitopes_to_fasta(all_epitopes, output_file)

# Example usage
if __name__ == "__main__":
    # Check if the CSV file path is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python batch_generate_epitopes.py [path to your csv file]")
        sys.exit(1)

    # Get the CSV file path from the command-line arguments
    csv_file_path = sys.argv[1]
    process_csv_file(csv_file_path)