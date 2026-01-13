def calculate_base_ratio(fasta_file):
    base_counts = {}
    total_bases = 0

    with open(fasta_file, "r") as file:
        lines = file.readlines()
        sequence = ""

        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    for base in sequence:
                        base = base.upper()  # 将碱基转换为大写
                        if base in base_counts:
                            base_counts[base] += 1
                        else:
                            base_counts[base] = 1
                        total_bases += 1
                    sequence = ""
            else:
                sequence += line

        for base, count in base_counts.items():
            ratio = count / total_bases
            print(f"{base}: {ratio:.6f}")

fasta_file = "/mnt/NFS/fengch/H99RNA.fna"
calculate_base_ratio(fasta_file)
