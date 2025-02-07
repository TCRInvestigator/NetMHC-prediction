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
    print(f"DEBUG: [fetch_latest_version] Response: {response[:200]}...")  # first 200 characters
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
    """
    Fetch the GenBank record for the given transcript ID and extract the CDS translation.
    """
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
    """
    Translate the mRNA sequence from the given start index until a stop codon.
    """
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
    """
    Determine the actual start codon by comparing translation from candidate 'ATG' positions
    to the expected protein sequence fetched from the transcript's GenBank record.
    
    Returns the trimmed mRNA sequence (starting at the actual start codon) and the candidate index.
    """
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
    """
    Translate the trimmed mRNA sequence into a protein until the first stop codon.
    """
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
    """
    For each versioned transcript, fetch the CDS protein from NIH and check if the
    amino acid (or protein length for stoploss mutations) at the mutation position matches
    the expected one extracted from the mutation info.
    """
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
    """
    Extract the WT amino acid and its 1-based mutation position from the protein change annotation.
    For example:
      - "p.L330Rfs*4" returns ("L", 330)
      - "p.G18_D19insGGGGGGG" returns ("G", 18)
      - "p.*1275delinsAGSWGWPGTTGCQPRARSLEGSWGNPSLLLDVCVTSVSPVLRWGISRAVVGQS*" returns ("*", 1275)
    """
    m = re.search(r'p\.([*A-Z])(\d+)', protein_change)
    if not m:
        raise ValueError(f"Unable to parse protein_change: {protein_change}")
    original_aa = m.group(1)
    mutation_position = int(m.group(2))  # 1-based position
    return original_aa, mutation_position

# --------------------------------------------------
# Epitope generation and filtering functions
# --------------------------------------------------

def generate_epitopes(protein_sequence, epitope_length=9):
    """Generate a list of overlapping peptides (epitopes) of the given length from a protein sequence."""
    epitopes = []
    for i in range(len(protein_sequence) - epitope_length + 1):
        epitopes.append(protein_sequence[i:i+epitope_length])
    return epitopes

def filter_unique_epitopes(wt_epitopes, mutated_epitopes):
    """Return only those peptides that appear in the mutated protein but not in the WT."""
    wt_set = set(wt_epitopes)
    unique = [ep for ep in mutated_epitopes if ep not in wt_set]
    return unique

def write_epitopes_to_fasta(epitope_entries, output_file):
    """Write the list of epitope entries to a FASTA file."""
    with open(output_file, 'w') as f:
        for entry in epitope_entries:
            f.write(entry + "\n")
    print(f"Epitopes have been written to {output_file}")

# --------------------------------------------------
# Mutation function (updated to assume trimmed mRNA)
# --------------------------------------------------

def apply_indel_mutation(trimmed_sequence, mutation_info):
    """
    Apply an indel mutation (deletion, insertion, or duplication) to the already trimmed mRNA sequence.
    The mutation_info string is in the format:
        gene:transcript:exon:dna_change:protein_change
    """
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
                mutation_position = mutation_position_1based - 1  # Convert to 0-based indexing
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

            # Trim the mRNA sequence to the actual start codon.
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

            wt_epitopes = generate_epitopes(wt_protein, 9)
            mutated_epitopes = generate_epitopes(mutated_protein, 9)
            unique_epitopes = filter_unique_epitopes(wt_epitopes, mutated_epitopes)

            for j, ep in enumerate(unique_epitopes, start=1):
                header = f"{mutation_info}_epitope_{j}"
                entry = f">{header}\n{ep}"
                all_epitope_entries.append(entry)

    base_name = os.path.splitext(os.path.basename(csv_file_path))[0]
    output_file = f"{base_name}_unique_epitopes.fasta"
    write_epitopes_to_fasta(all_epitope_entries, output_file)

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
