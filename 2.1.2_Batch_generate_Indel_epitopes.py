import pandas as pd
import subprocess
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
import os
import sys
import re

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

# --------------------------------------------------
# Basic functions (unchanged from the original script)
# --------------------------------------------------

def clean_dna_sequence(data):
    """Remove header lines and whitespace from a FASTA-formatted sequence."""
    lines = data.strip().split('\n')
    sequence_lines = [line for line in lines if not line.startswith('>')]
    return ''.join(sequence_lines).replace(' ', '').replace('\n', '')

def fetch_latest_version(transcript_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={transcript_id}"
    response = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    match = re.search(r"<Id>(\d+)<\/Id>", response)
    if match:
        versioned_id = match.group(1)
    else:
        raise ValueError(f"Could not determine the ID for transcript {transcript_id}")

    summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={versioned_id}"
    summary_response = subprocess.run(['curl', '-s', summary_url], capture_output=True, text=True, check=True).stdout

    version_match = re.search(r"gi\|\d+\|ref\|(NM_\d+\.\d+)\|", summary_response)
    if version_match:
        versioned_transcript_id = version_match.group(1)
        return versioned_transcript_id
    else:
        raise ValueError(f"Could not determine the versioned transcript ID for ID {versioned_id}")

def fetch_mrna_sequence(transcript_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={transcript_id}&rettype=fasta&retmode=text"
    nucleotide_record = subprocess.run(['curl', '-s', url], capture_output=True, text=True, check=True).stdout
    return clean_dna_sequence(nucleotide_record)

def trim_to_first_atg(sequence):
    atg_index = sequence.find('ATG')
    if atg_index == -1:
        raise ValueError("No ATG start codon found in the sequence.")
    return sequence[atg_index:], atg_index

def translate_mrna_to_protein(sequence):
    """
    Translate the mRNA sequence into a protein.
    The mRNA is first trimmed to the first ATG; then translation
    proceeds codon-by-codon until the first stop codon.
    """
    sequence, _ = trim_to_first_atg(sequence)
    protein_sequence = ""
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) < 3:
            break
        amino_acid = str(Seq(codon).translate())
        if amino_acid == '*':  # Stop codon encountered
            break
        protein_sequence += amino_acid
    return protein_sequence

def get_all_versioned_transcripts(transcript_id, latest_version):
    version_list = []
    version_number = int(latest_version.split('.')[-1])
    for version in range(1, version_number + 1):
        versioned_transcript = f"{transcript_id}.{version}"
        version_list.append(versioned_transcript)
    return version_list

def find_matching_version(transcript_id, mutation_position, original_aa):
    """
    Loop over versioned transcripts until one is found where the
    WT protein (translated from the mRNA trimmed to the first ATG)
    has the expected WT amino acid at the given (0-based) mutation position.
    
    For stoploss mutations (where original_aa == "*"), we require that
    the wild-type protein length equals the mutation position (i.e. the stop codon is absent).
    """
    latest_version = fetch_latest_version(transcript_id)
    versioned_transcripts = get_all_versioned_transcripts(transcript_id, latest_version)

    for versioned_transcript in versioned_transcripts:
        try:
            mrna_sequence = fetch_mrna_sequence(versioned_transcript)
            protein_sequence = translate_mrna_to_protein(mrna_sequence)
            if original_aa == "*":
                # For stoploss, the WT protein (which does not include the stop codon)
                # should have length exactly equal to the 0-based mutation position.
                if len(protein_sequence) == mutation_position:
                    return versioned_transcript, mrna_sequence, protein_sequence
            else:
                if len(protein_sequence) > mutation_position and protein_sequence[mutation_position] == original_aa:
                    return versioned_transcript, mrna_sequence, protein_sequence
        except ValueError as e:
            print(f"Error processing {versioned_transcript}: {e}")
            continue
    raise ValueError(f"No matching version found for transcript ID {transcript_id} with specified mutation.")

def apply_indel_mutation(sequence, mutation_info):
    """
    Apply an indel mutation (deletion, insertion, or duplication) to the mRNA sequence.
    The input mRNA sequence is trimmed to the first ATG.
    The mutation_info string is expected to have the format:
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

    trimmed_sequence, atg_index = trim_to_first_atg(sequence)
    adjusted_mutation_position = mutation_position
    if adjusted_mutation_position < 0:
        raise ValueError(f"Mutation position {mutation_position} is before the first ATG in the sequence.")

    if mutation_type == 'del':
        mutated_sequence = trimmed_sequence[:adjusted_mutation_position] + trimmed_sequence[adjusted_mutation_position + delete_length:]
    elif mutation_type in ['ins', 'dup']:
        mutated_sequence = trimmed_sequence[:adjusted_mutation_position + 1] + inserted_nucleotides + trimmed_sequence[adjusted_mutation_position + 1:]
    else:
        raise ValueError(f"Unsupported mutation type: {mutation_type}")
    return mutated_sequence

# --------------------------------------------------
# New helper function for mutation info extraction
# --------------------------------------------------

def extract_mutation_info(protein_change):
    """
    Extract the WT amino acid and its position from the protein change annotation.
    Examples:
      - For "p.L330Rfs*4", returns ("L", 330)
      - For "p.G18_D19insGGGGGGG", returns ("G", 18)
      - For "p.*1275delinsAGSWGWPGTTGCQPRARSLEGSWGNPSLLLDVCVTSVSPVLRWGISRAVVGQS*", returns ("*", 1275)
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

def write_epitopes_to_fasta(epitopes, output_file):
    """Write the list of epitopes to a FASTA file."""
    with open(output_file, 'w') as f:
        for i, ep in enumerate(epitopes):
            f.write(f">epitope_{i+1}\n{ep}\n")
    print(f"Epitopes have been written to {output_file}")

# --------------------------------------------------
# Main processing function
# --------------------------------------------------

def process_csv_file(csv_file_path):
    df = pd.read_csv(csv_file_path)
    all_unique_epitopes = []
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

            # Extract the WT amino acid and its 1-based mutation position.
            try:
                original_aa, mutation_position_1based = extract_mutation_info(protein_change)
                mutation_position = mutation_position_1based - 1  # Convert to 0-based
            except Exception as e:
                print(f"Error extracting WT info from {protein_change}: {e}")
                unmatched_mutations.append(mutation_info)
                continue

            try:
                matching_transcript, mrna_sequence, wt_protein = find_matching_version(transcript, mutation_position, original_aa)
            except ValueError as e:
                print(f"Error finding matching version for {mutation_info}: {e}")
                unmatched_mutations.append(mutation_info)
                continue

            try:
                mutated_mrna = apply_indel_mutation(mrna_sequence, mutation_info)
            except ValueError as e:
                print(f"Error applying indel mutation for {mutation_info}: {e}")
                unmatched_mutations.append(mutation_info)
                continue

            # Translate both WT and mutated mRNA sequences (both trimmed to first ATG).
            mutated_protein = translate_mrna_to_protein(mutated_mrna)
            # wt_protein is already obtained

            wt_epitopes = generate_epitopes(wt_protein, 9)
            mutated_epitopes = generate_epitopes(mutated_protein, 9)

            unique_epitopes = filter_unique_epitopes(wt_epitopes, mutated_epitopes)

            for ep in unique_epitopes:
                header = f"{mutation_info}_unique_epitope"
                all_unique_epitopes.append(f">{header}\n{ep}")

    base_name = os.path.splitext(os.path.basename(csv_file_path))[0]
    output_file = f"{base_name}_unique_epitopes.fasta"
    write_epitopes_to_fasta(all_unique_epitopes, output_file)

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
