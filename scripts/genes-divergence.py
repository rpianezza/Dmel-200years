# Takes as input the folder where the BUSCO subfolders are stored (output of "busco.sh").
# Outputs a fasta file for each BUSCO entry shared by all the species and a list of the shared BUSCO entries 

import os
import argparse
import subprocess

def read_entries(file_path):
    entries = set()
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if line.startswith("#"):
                pass
            elif columns[1] == 'Complete':
                entries.add(columns[0])
    return entries

def find_common_entries(root_folder, species1, species2):
    common_entries = set()

    for species in [species1, species2]:
        subfolder_path = os.path.join(root_folder, species)
        subsubfolder_path = os.path.join(subfolder_path, 'run_diptera_odb10')

        if os.path.exists(subsubfolder_path):
            current_entries = read_entries(os.path.join(subsubfolder_path, 'full_table.tsv'))
            if not common_entries:
                common_entries.update(current_entries)
            else:
                common_entries.intersection_update(current_entries)

    return list(common_entries)

def extract_sequences(gene, root_folder, species1, species2, output):
    for species in [species1, species2]:
        subfolder_path = os.path.join(root_folder, species)
        subsubfolder_path = os.path.join(subfolder_path, 'run_diptera_odb10')
        if os.path.exists(subsubfolder_path):
            fasta_file_path = os.path.join(subsubfolder_path, 'busco_sequences', 'single_copy_busco_sequences', f'{gene}.fna')
            if os.path.exists(fasta_file_path):
                with open(fasta_file_path, 'r') as fasta_file, open(output, 'a') as output_file:
                    for record in fasta_file:
                        if record.startswith('>'):
                            output_file.write('>'+species+'\n')
                        else:
                            output_file.write(record)

def run_muscle(input_folder, output_folder):
    # Iterate over files in the input folder
    for filename in os.listdir(input_folder):
        input_file_path = os.path.join(input_folder, filename)
        
        # Check if the item in the folder is a file
        if os.path.isfile(input_file_path):
            output_file_path = os.path.join(output_folder, filename + ".MSA")
            
            # Run MUSCLE on the input file
            subprocess.run(["muscle", "-in", input_file_path, "-out", output_file_path])

def calculate_divergence(msa_file):
    with open(msa_file, 'r') as file:
        lines = file.readlines()
        seq1 = lines[1].strip()  # Assuming the sequences are on the 2nd and 4th lines
        seq2 = lines[3].strip()

    num_differences = sum(a != b for a, b in zip(seq1, seq2))
    alignment_length = len(seq1)

    percent_divergence = (num_differences / alignment_length) * 100
    return percent_divergence

def defragment(fasta, out):
    with open(fasta, 'r') as fasta_file, open(out, 'w') as output_file:
        for i, line in enumerate(fasta_file):
            line = line.strip()
            if line.startswith(">"):
                if i>0:
                    output_file.write("\n")
                output_file.write(line + "\n")
            else:
                output_file.write(line)

def main():
    parser = argparse.ArgumentParser(description='Extract and find common entries from TSV files.')
    parser.add_argument('speciesA', help='Name of species 1')
    parser.add_argument('speciesB', help='Name of species 2')
    parser.add_argument('input_folder', help='Path to the input folder containing subfolders with TSV files')
    parser.add_argument('output', help='Path to the output folder')

    args = parser.parse_args()

    # Create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Create subfolders "fasta" and "muscle" inside the output folder
    fasta_folder = os.path.join(args.output, 'fasta')
    muscle_folder = os.path.join(args.output, 'muscle')
    for folder in [fasta_folder, muscle_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    result = find_common_entries(args.input_folder, args.speciesA, args.speciesB)

    for gene in result:
        print(gene)
        extract_sequences(gene, args.input_folder, args.speciesA, args.speciesB, fasta_folder+'/'+gene+'.fasta')

    run_muscle(fasta_folder, muscle_folder)

    output_filename = os.path.join(args.output, f"divergence-{args.speciesA}-{args.speciesB}.tsv")
    with open(output_filename, 'w') as output:
        output.write("gene\tdivergence\n")

        for filename in os.listdir(muscle_folder):
            msa = os.path.join(muscle_folder, filename)
            defragment(msa, msa+".defrag")
        
        for filename in os.listdir(muscle_folder):
            if filename.endswith(".defrag"):  # Assuming all MSA files have the .MSA extension
                msa_file = os.path.join(muscle_folder, filename)
                divergence = calculate_divergence(msa_file)
                        
                # Extract gene name from filename
                gene = filename.replace(".fasta.MSA.defrag", "")
                        
                # Write gene name and divergence to output file
                output.write(f"{gene}\t{divergence}\n")

if __name__ == "__main__":
    main()