# Practica 1
# Equipo: Los bipedos
# Integrantes:
#  - Francisco Contreras Ibarra
#  - Jose Ethan Ortega Gonzalez
#  - Alvaro Ramirez Lopez

import sys

def translate_codon(codon):
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    return codon_table.get(codon, 'X')

def dna_to_protein(dna_sequence):
    protein_sequence = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        amino_acid = translate_codon(codon)
        if amino_acid != '*':
            protein_sequence += amino_acid
        else:
            break
    return protein_sequence

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python p1_losBipedos.py <archivo_entrada>.fasta <archivo_salida>.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, "r") as f_in:
        dna_sequence = "".join(line.strip() for line in f_in.readlines()[1:])

    protein_sequence = dna_to_protein(dna_sequence)

    with open(output_file, "w") as f_out:
        f_out.write(f">{input_file}\n")
        f_out.write(protein_sequence)
