def calculate_aa_composition(fasta_file):
    aa_counts = {}
    total_aa_count = 0

    with open(fasta_file, 'r') as file:
        current_sequence = ''

        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    for aa in current_sequence:
                        aa_counts[aa] = aa_counts.get(aa, 0) + 1
                        total_aa_count += 1
                current_sequence = ''
            else:
                current_sequence += line

        for aa in current_sequence:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
            total_aa_count += 1

    aa_composition = {}
    for aa, count in aa_counts.items():
        composition = count / total_aa_count
        aa_composition[aa] = round(composition, 6)  # Keep 6 decimal places

    return aa_composition


fasta_file = '/mnt/NFS/fengch/H99protein.faa'
result = calculate_aa_composition(fasta_file)

for aa, composition in result.items():
    print(f"Amino Acid: {aa}, Composition: {composition:.6f}")  # Print with 6 decimal plac
