import pandas as pd
import subprocess
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
import os
import sys
import re
from io import StringIO
from Bio import SeqIO
import requests  # For UniProt queries

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

# --------------------------------------------------
# Global API Key
# --------------------------------------------------
API_KEY = "60e3f2a2a13429eb6156a6a1989441e29d08"

# --------------------------------------------------
# Basic functions
# --------------------------------------------------

def clean_dna_sequence(data):
    """Remove header lines and whitespace from a FASTA-formatted sequence."""
    lines = data.strip().split('\n')
    sequence_lines = [line for line in lines if not line.startswith('>')]
    return ''.join(sequence_lines).replace(' ', '').replace('\n', '')

def fetch_latest_version(transcript_id):
    print(f"DEBUG: [fetch_latest_version] Transcript ID: {transcript_id}")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={transcript_id}&api_key={API_KEY}"
    print(f"DEBUG: [fetch_latest_version] URL: {url}")
    try:
        response = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError as e:
        print(f"DEBUG: [fetch_latest_version] Curl command failed when accessing URL: {url}")
        print(f"DEBUG: Error: {e}")
        raise
    print(f"DEBUG: [fetch_latest_version] Response: {response[:200]}...")
    match = re.search(r"<Id>(\d+)<\/Id>", response)
    if match:
        versioned_id = match.group(1)
        print(f"DEBUG: [fetch_latest_version] Found versioned ID: {versioned_id}")
    else:
        raise ValueError(f"Could not determine the ID for transcript {transcript_id}")
    
    summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={versioned_id}&api_key={API_KEY}"
    print(f"DEBUG: [fetch_latest_version] Summary URL: {summary_url}")
    try:
        summary_response = subprocess.run(['curl', '-s', summary_url], capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError as e:
        print(f"DEBUG: [fetch_latest_version] Curl command failed when accessing URL: {summary_url}")
        print(f"DEBUG: Error: {e}")
        raise
    print(f"DEBUG: [fetch_latest_version] Summary Response: {summary_response[:200]}...")
    version_match = re.search(r"gi\|\d+\|ref\|(NM_\d+\.\d+)\|", summary_response)
    if version_match:
        versioned_transcript_id = version_match.group(1)
        print(f"DEBUG: [fetch_latest_version] Versioned transcript ID: {versioned_transcript_id}")
        return versioned_transcript_id
    else:
        raise ValueError(f"Could not determine the versioned transcript ID for ID {versioned_id}")

def fetch_mrna_sequence(transcript_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={transcript_id}&rettype=fasta&retmode=text&api_key={API_KEY}"
    nucleotide_record = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    return clean_dna_sequence(nucleotide_record)

# --------------------------------------------------
# Helper functions for start codon determination
# --------------------------------------------------

def fetch_protein_sequence_from_transcript(transcript_id):
    print(f"DEBUG: [fetch_protein_sequence_from_transcript] Fetching GenBank record for transcript: {transcript_id}")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={transcript_id}&rettype=gb&retmode=text&api_key={API_KEY}"
    genbank_text = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    record = SeqIO.read(StringIO(genbank_text), "genbank")
    print(f"DEBUG: [fetch_protein_sequence_from_transcript] Found GenBank record with {len(record.features)} features.")
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            translation = feature.qualifiers["translation"][0]
            print(f"DEBUG: [fetch_protein_sequence_from_transcript] Found CDS with translation length {len(translation)}")
            return translation
    raise ValueError(f"No CDS translation found in GenBank record for {transcript_id}")

def translate_from_start(sequence, start_index):
    protein = ""
    for i in range(start_index, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) < 3:
            break
        aa = str(Seq(codon).translate())
        if aa == "*":
            break
        protein += aa
    print(f"DEBUG: [translate_from_start] Translated protein from index {start_index}: {protein}")
    return protein

def trim_to_first_start_codon(mrna_sequence, transcript_id):
    print(f"DEBUG: [trim_to_first_start_codon] Received mRNA sequence of length {len(mrna_sequence)} for transcript {transcript_id}")
    expected_protein = fetch_protein_sequence_from_transcript(transcript_id)
    print(f"DEBUG: [trim_to_first_start_codon] Expected protein (length {len(expected_protein)}): {expected_protein}")
    pos = 0
    while True:
        candidate = mrna_sequence.find("ATG", pos)
        if candidate == -1:
            print("DEBUG: [trim_to_first_start_codon] No further 'ATG' found in the mRNA sequence.")
            break
        candidate_protein = translate_from_start(mrna_sequence, candidate)
        print(f"DEBUG: [trim_to_first_start_codon] Trying ATG at position {candidate}: candidate protein (length {len(candidate_protein)})")
        if candidate_protein == expected_protein:
            print(f"DEBUG: [trim_to_first_start_codon] Found matching start codon at position {candidate}")
            return mrna_sequence[candidate:], candidate
        pos = candidate + 1
    raise ValueError("No start codon found that produces the expected protein sequence.")

# --------------------------------------------------
# Translation function (assumes input is already trimmed)
# --------------------------------------------------

def translate_mrna_to_protein(trimmed_sequence):
    protein_sequence = ""
    for i in range(0, len(trimmed_sequence) - 2, 3):
        codon = trimmed_sequence[i:i+3]
        if len(codon) < 3:
            break
        aa = str(Seq(codon).translate())
        if aa == "*":
            break
        protein_sequence += aa
    print(f"DEBUG: [translate_mrna_to_protein] Translated protein: {protein_sequence}")
    return protein_sequence

# --------------------------------------------------
# Version helper
# --------------------------------------------------

def get_all_versioned_transcripts(transcript_id, latest_version):
    print(f"DEBUG: [get_all_versioned_transcripts] Transcript ID: {transcript_id}, Latest version: {latest_version}")
    version_list = []
    try:
        version_number = int(latest_version.split('.')[-1])
    except Exception as e:
        raise ValueError(f"Error parsing latest version number from {latest_version}: {e}")
    for version in range(1, version_number + 1):
        versioned_transcript = f"{transcript_id}.{version}"
        version_list.append(versioned_transcript)
    print(f"DEBUG: [get_all_versioned_transcripts] Constructed version list: {version_list}")
    return version_list

# --------------------------------------------------
# Matching function (updated to use NIH protein sequence)
# --------------------------------------------------

def find_matching_version(transcript_id, mutation_position, original_aa, protein_change):
    print(f"DEBUG: [find_matching_version] Starting search for transcript {transcript_id}")
    latest_version = fetch_latest_version(transcript_id)
    versioned_transcripts = get_all_versioned_transcripts(transcript_id, latest_version)
    for versioned_transcript in versioned_transcripts:
        print(f"DEBUG: [find_matching_version] Trying version: {versioned_transcript}")
        try:
            mrna_sequence = fetch_mrna_sequence(versioned_transcript)
            protein_sequence = fetch_protein_sequence_from_transcript(versioned_transcript)
            print(f"DEBUG: [find_matching_version] For version {versioned_transcript}, protein length: {len(protein_sequence)}")
            if original_aa == "*":
                if len(protein_sequence) == mutation_position:
                    print(f"DEBUG: [find_matching_version] Stoploss match found for {versioned_transcript}")
                    return versioned_transcript, mrna_sequence, protein_sequence
                else:
                    print(f"DEBUG: [find_matching_version] Stoploss match failed for {versioned_transcript}: protein length {len(protein_sequence)} vs mutation position {mutation_position}")
            else:
                if len(protein_sequence) > mutation_position:
                    print(f"DEBUG: [find_matching_version] Checking residue at position {mutation_position} in {versioned_transcript}: expected {original_aa}, found {protein_sequence[mutation_position]}")
                    if protein_sequence[mutation_position] == original_aa:
                        print(f"DEBUG: [find_matching_version] Match found for {versioned_transcript}")
                        return versioned_transcript, mrna_sequence, protein_sequence
                    else:
                        print(f"DEBUG: [find_matching_version] Residue mismatch for {versioned_transcript}")
                else:
                    print(f"DEBUG: [find_matching_version] Protein sequence too short for {versioned_transcript}")
        except ValueError as e:
            print(f"Error processing {versioned_transcript}: {e}")
            continue
    raise ValueError(f"No matching version found for transcript ID {transcript_id} with specified mutation.")

# --------------------------------------------------
# Helper function for mutation info extraction
# --------------------------------------------------

def extract_mutation_info(protein_change):
    m = re.search(r'p\.([*A-Z])(\d+)', protein_change)
    if not m:
        raise ValueError(f"Unable to parse protein_change: {protein_change}")
    original_aa = m.group(1)
    mutation_position = int(m.group(2))
    return original_aa, mutation_position

# --------------------------------------------------
# Epitope generation and filtering functions
# --------------------------------------------------

# Function to check if a peptide exists in the UniProt human proteome
def is_in_human_proteome(peptide):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={peptide}+AND+reviewed:true+AND+organism:9606"
    response = requests.get(url)
    return response.status_code == 200 and 'accession' in response.text

# Generate 15-mer peptides
def generate_epitopes(protein_sequence, epitope_length=15):
    return [protein_sequence[i:i+epitope_length] for i in range(len(protein_sequence) - epitope_length + 1)]

def filter_human_proteome(epitopes):
    query_size = 100
    epitopes = list(set(epitopes))
    human_proteome_hits = set()
    for i in range(0, len(epitopes), query_size):
        batch = epitopes[i:i+query_size]
        query = " OR ".join(batch)
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}+AND+reviewed:true+AND+organism:9606"
        response = requests.get(url)
        if response.status_code == 200:
            for epitope in batch:
                if epitope in response.text:
                    human_proteome_hits.add(epitope)
    return human_proteome_hits

def filter_unique_epitopes(wt_epitopes, mutated_epitopes):
    wt_set = set(wt_epitopes)
    unique = [ep for ep in mutated_epitopes if ep not in wt_set]
    return unique

def write_epitopes_to_fasta(epitope_entries, output_file):
    with open(output_file, 'w') as f:
        for entry in epitope_entries:
            f.write(entry + "\n")
    print(f"Epitopes have been written to {output_file}")

# --------------------------------------------------
# Mutation function (updated to assume trimmed mRNA)
# --------------------------------------------------

def apply_indel_mutation(trimmed_sequence, mutation_info):
    gene, transcript, exon, dna_change, protein_change = mutation_info.split(':')
    if 'del' in dna_change:
        mutation_type = 'del'
        if '_' in dna_change:
            start_pos, end_pos = map(int, re.findall(r'\d+', dna_change.split('del')[0]))
            mutation_position = start_pos - 1
            delete_length = end_pos - start_pos + 1
        else:
            mutation_position = int(''.join(filter(str.isdigit, dna_change.split('del')[0]))) - 1
            delete_length = 1
    elif 'ins' in dna_change:
        mutation_type = 'ins'
        start_pos, end_pos = map(int, re.findall(r'\d+', dna_change.split('ins')[0]))
        mutation_position = start_pos - 1
        inserted_nucleotides = dna_change.split('ins')[-1]
    elif 'dup' in dna_change:
        mutation_type = 'dup'
        mutation_position = int(''.join(filter(str.isdigit, dna_change.split('dup')[0]))) - 1
        inserted_nucleotides = dna_change.split('dup')[-1]
    else:
        raise ValueError(f"Unsupported mutation type in {dna_change}")
    if mutation_position < 0:
        raise ValueError(f"Mutation position {mutation_position} is before the start codon in the trimmed sequence.")
    if mutation_type == 'del':
        mutated_sequence = trimmed_sequence[:mutation_position] + trimmed_sequence[mutation_position + delete_length:]
    elif mutation_type in ['ins', 'dup']:
        mutated_sequence = trimmed_sequence[:mutation_position + 1] + inserted_nucleotides + trimmed_sequence[mutation_position + 1:]
    else:
        raise ValueError(f"Unsupported mutation type: {mutation_type}")
    print(f"DEBUG: [apply_indel_mutation] Applied {mutation_type} mutation at position {mutation_position}.")
    return mutated_sequence

# --------------------------------------------------
# Main processing function
# --------------------------------------------------

def process_csv_file(csv_file_path):
    df = pd.read_csv(csv_file_path)
    all_epitope_entries = []
    all_epitopes = []  # List to store all generated epitopes
    mutation_map = []  # List to store corresponding mutation info for each epitope
    unmatched_mutations = []

    for index, row in df.iterrows():
        if pd.isna(row['AAChange']):
            continue
        aa_changes = row['AAChange'].split(',')
        for mutation_info in aa_changes:
            try:
                gene, transcript, exon, dna_change, protein_change = mutation_info.split(':')
            except Exception as e:
                print(f"Error parsing mutation info {mutation_info}: {e}")
                unmatched_mutations.append(mutation_info)
                continue
            try:
                original_aa, mutation_position_1based = extract_mutation_info(protein_change)
                mutation_position = mutation_position_1based - 1
            except Exception as e:
                print(f"Error extracting WT info from {protein_change}: {e}")
                unmatched_mutations.append(mutation_info)
                continue
            try:
                matching_transcript, mrna_sequence, wt_protein = find_matching_version(
                    transcript, mutation_position, original_aa, protein_change
                )
            except ValueError as e:
                print(f"Error finding matching version for {mutation_info}: {e}")
                unmatched_mutations.append(mutation_info)
                continue
            try:
                trimmed_mrna, start_index = trim_to_first_start_codon(mrna_sequence, matching_transcript)
            except ValueError as e:
                print(f"Error trimming mRNA for {mutation_info}: {e}")
                unmatched_mutations.append(mutation_info)
                continue
            try:
                mutated_mrna = apply_indel_mutation(trimmed_mrna, mutation_info)
            except ValueError as e:
                print(f"Error applying indel mutation for {mutation_info}: {e}")
                unmatched_mutations.append(mutation_info)
                continue
            try:
                mutated_protein = translate_mrna_to_protein(mutated_mrna)
            except ValueError as e:
                print(f"Error translating mutated mRNA for {mutation_info}: {e}")
                unmatched_mutations.append(mutation_info)
                continue
            wt_epitopes = generate_epitopes(wt_protein, 15)
            mutated_epitopes = generate_epitopes(mutated_protein, 15)
            unique_epitopes = filter_unique_epitopes(wt_epitopes, mutated_epitopes)
            for ep in unique_epitopes:
                all_epitopes.append(ep)
                mutation_map.append(mutation_info)
    # Bulk filtering via UniProt
    print("Querying UniProt for filtering...")
    human_proteome_hits = filter_human_proteome(all_epitopes)
    
    # Separate into kept and removed pairs
    kept_pairs = []
    removed_pairs = []
    for m_info, ep in zip(mutation_map, all_epitopes):
        if ep in human_proteome_hits:
            removed_pairs.append((m_info, ep))
        else:
            kept_pairs.append((m_info, ep))
    
    # Generate FASTA entries for kept epitopes using mutation-specific numbering
    epitope_counter = {}
    final_epitope_entries = []
    for m_info, ep in kept_pairs:
        if m_info not in epitope_counter:
            epitope_counter[m_info] = 1
        else:
            epitope_counter[m_info] += 1
        header = f"{m_info}_epitope_{epitope_counter[m_info]}"
        entry = f">{header}\n{ep}"
        final_epitope_entries.append(entry)
    
    # Generate FASTA entries for removed epitopes using mutation-specific numbering
    removed_epitope_counter = {}
    removed_epitope_entries = []
    for m_info, ep in removed_pairs:
        if m_info not in removed_epitope_counter:
            removed_epitope_counter[m_info] = 1
        else:
            removed_epitope_counter[m_info] += 1
        header = f"{m_info}_removed_epitope_{removed_epitope_counter[m_info]}"
        entry = f">{header}\n{ep}"
        removed_epitope_entries.append(entry)
    
    base_name = os.path.splitext(os.path.basename(csv_file_path))[0]
    kept_output_file = f"{base_name}_unique_epitopes.fasta"
    removed_output_file = f"{base_name}_removed_epitopes.fasta"
    
    write_epitopes_to_fasta(final_epitope_entries, kept_output_file)
    write_epitopes_to_fasta(removed_epitope_entries, removed_output_file)
    
    unmatched_output_file = f"{base_name}_unmatched_mutations.csv"
    pd.DataFrame(unmatched_mutations, columns=['UnmatchedMutation']).to_csv(unmatched_output_file, index=False)
    print(f"Unmatched mutations have been written to {unmatched_output_file}")

# --------------------------------------------------
# Main execution
# --------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python batch_generate_epitopes.py [path to your csv file]")
        sys.exit(1)
    csv_file_path = sys.argv[1]
    process_csv_file(csv_file_path)
