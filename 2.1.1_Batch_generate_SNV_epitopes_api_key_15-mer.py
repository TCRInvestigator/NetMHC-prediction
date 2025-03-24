import pandas as pd
import subprocess
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
import os
import sys
import subprocess
import re

# Function to clean the DNA sequence by removing non-sequence lines
def clean_dna_sequence(data):
    lines = data.strip().split('\n')
    # Remove any header or annotation lines (those starting with '>')
    sequence_lines = [line for line in lines if not line.startswith('>')]
    # Join the remaining lines to form the clean sequence
    return ''.join(sequence_lines).replace(' ', '').replace('\n', '')

# Function to fetch the latest version number of the transcript ID
def fetch_latest_version(transcript_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={transcript_id}"
    response = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    match = re.search(r"<Id>(\d+)<\/Id>", response)
    if match:
        versioned_id = match.group(1)
    else:
        raise ValueError(f"Could not determine the ID for transcript {transcript_id}")

    # Step 2: Fetch the summary for that ID to get the versioned transcript
    summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={versioned_id}"
    summary_response = subprocess.run(['curl', '-s', summary_url], capture_output=True, text=True, check=True).stdout

    # Extract versioned accession from the summary (from the "Extra" line)
    version_match = re.search(r"gi\|\d+\|ref\|(NM_\d+\.\d+)\|", summary_response)
    if version_match:
        versioned_transcript_id = version_match.group(1)
        return versioned_transcript_id
    else:
        raise ValueError(f"Could not determine the versioned transcript ID for ID {versioned_id}")

# Function to fetch protein sequence directly from NCBI
def fetch_protein_sequence(transcript_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={transcript_id}&rettype=gb&retmode=text"
    nucleotide_record = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    match = re.search(r'/translation="([A-Z\s]+)"', nucleotide_record)
    if match:
        protein_sequence = "".join(match.group(1).split())
        return protein_sequence
    else:
        raise ValueError(f"Could not fetch protein sequence for transcript ID {transcript_id}")

# Function to get all versioned transcript IDs
def get_all_versioned_transcripts(transcript_id, latest_version):
    version_list = []
    version_number = int(latest_version.split('.')[-1])
    for version in range(1, version_number + 1):
        versioned_transcript = f"{transcript_id}.{version}"
        version_list.append(versioned_transcript)
    return version_list

# Function to validate which version matches the mutation information
def find_matching_version(transcript_id, mutation_position, original_aa):
    latest_version = fetch_latest_version(transcript_id)
    versioned_transcripts = get_all_versioned_transcripts(transcript_id, latest_version)

    for versioned_transcript in versioned_transcripts:
        try:
            protein_sequence = fetch_protein_sequence(versioned_transcript)
            if len(protein_sequence) > mutation_position and protein_sequence[mutation_position] == original_aa:
                return versioned_transcript, protein_sequence
        except ValueError:
            continue
    raise ValueError(f"No matching version found for transcript ID {transcript_id} with specified mutation.")

# Function to apply mutation and generate epitopes
def generate_epitopes(protein_sequence, mutation_info):
    # Parse mutation information
    gene, transcript, exon, dna_change, protein_change = mutation_info.split(':')
    
    # Extract the position number from protein_change, e.g., extract 171 from p.Q171R
    mutation_position_aa = int(''.join(filter(str.isdigit, protein_change))) - 1  # Convert to 0-based index
    new_aa = protein_change[-1]

    # Apply the mutation directly to the protein sequence
    mutated_protein_sequence = (protein_sequence[:mutation_position_aa] + new_aa + protein_sequence[mutation_position_aa + 1:])

    protein_len = len(mutated_protein_sequence)

    wildtype_epitopes = []
    mutated_epitopes = []
    epitope_length = 15  # Length of the epitope

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

        wildtype_epitope_sequence = protein_sequence[start:end]
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
    unmatched_mutations = []

    # Process each mutation and generate epitopes
    for index, row in df.iterrows():
        aa_changes = row['AAChange'].split(',')
        for mutation_info in aa_changes:
            gene, transcript, exon, dna_change, protein_change = mutation_info.split(':')
            
            # Extract the position number from protein_change, e.g., extract 171 from p.Q171R
            mutation_position_aa = int(''.join(filter(str.isdigit, protein_change))) - 1  # Convert to 0-based index
            original_aa = protein_change[2]  # Extract the original amino acid

            # Find the matching version of the transcript
            try:
                matching_transcript, protein_sequence = find_matching_version(transcript, mutation_position_aa, original_aa)
            except ValueError:
                unmatched_mutations.append(mutation_info)
                continue

            # Generate epitopes using the matching version
            wildtype_epitopes, mutated_epitopes = generate_epitopes(protein_sequence, mutation_info)
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

    # Write unmatched mutations to a CSV file
    unmatched_output_file = f"{base_name}_unmatched_mutations.csv"
    unmatched_df = pd.DataFrame(unmatched_mutations, columns=['UnmatchedMutation'])
    unmatched_df.to_csv(unmatched_output_file, index=False)
    print(f"Unmatched mutations have been written to {unmatched_output_file}")

# Example usage
if __name__ == "__main__":
    # Check if the CSV file path is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python batch_generate_epitopes.py [path to your csv file]")
        sys.exit(1)

    # Get the CSV file path from the command-line arguments
    csv_file_path = sys.argv[1]
    process_csv_file(csv_file_path)
