import random
from codon_dict import *
def dna_to_aa(dna_sequence):
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    if len(dna_sequence) % 3 != 0:
        raise ValueError("DNA sequence length must be a multiple of 3")
    amino_acids = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if codon in codon_table:
            amino_acids.append(codon_table[codon])
        else:
            raise ValueError("Unknown codon: {}".format(codon))
    
    return ''.join(amino_acids)

def optimize_gene_sequence(aa_sequence, s, selected_host_organism):
    if selected_host_organism == 'Phaseolus vulgaris':
        codon_bias_dict=codon_bias_dict_phaseolus_vulgaris
    if selected_host_organism == 'lens culinaris':
        codon_bias_dict=codon_bias_dict_lentil_culinaris
    if selected_host_organism == 'Cajanus cajan':
        codon_bias_dict=codon_bias_dict_cajanus_cajan
    if selected_host_organism == 'Vigna mungo':
        codon_bias_dict=codon_bias_dict_vigna_mungo
    if selected_host_organism == 'Vigna radiata':
        codon_bias_dict=codon_bias_dict_vigna_radiata
    if selected_host_organism == 'Lathyrus sativus':
        codon_bias_dict=codon_bias_dict_lathyrus_sativus
    if selected_host_organism == 'Pisum sativum':
        codon_bias_dict=codon_bias_dict_pisum_sativum
    if selected_host_organism == 'Glycine max':
        codon_bias_dict=codon_bias_dict_glycinemax
    if selected_host_organism == 'Cicer arientum':
        codon_bias_dict=codon_bias_dict_cicer_arientum
        

    optimized_dna_sequence = ""
    modify_seq_nospace=""
    for aa in aa_sequence:
        if aa in codon_bias_dict:
            codon_frequencies = codon_bias_dict[aa]
            # When s equals to 0, only the most frequently used codon is selected (full optimization)
            if s == 0:
                optimal_codon = max(codon_frequencies, key=codon_frequencies.get)
                optimized_dna_sequence += optimal_codon
                modify_seq_nospace += optimal_codon
            # When s equals to a large number, the least frequently used codons are selected (anti-optimization)
            elif s == 64:
                least_frequent_codon = min(codon_frequencies, key=codon_frequencies.get)
                optimized_dna_sequence += least_frequent_codon
                modify_seq_nospace += least_frequent_codon
            # When s equals to 1, each codon is equally weighted, resulting in random selection
            elif s == 1:
                selected_codon = random.choice(list(codon_frequencies.keys()))
                optimized_dna_sequence += selected_codon
                modify_seq_nospace += selected_codon
            else:
                upper_boundaries = [sum(list(codon_frequencies.values())[:i + 1]) ** s for i in range(len(codon_frequencies))]
                probabilities = [upper_boundaries[i] - upper_boundaries[i - 1] if i > 0 else upper_boundaries[i] for i in range(len(upper_boundaries))]
                selected_codon = random.choices(list(codon_frequencies.keys()), weights=probabilities, k=1)[0]
                    # optimized_dna_sequence += selected_codon + " "
                optimized_dna_sequence += selected_codon
                modify_seq_nospace += (selected_codon)
        else:
            optimized_dna_sequence += "*" # Add space if amino acid not found in codon bias (e.g., stop codon)
    return optimized_dna_sequence,modify_seq_nospace  # Strip trailing whitespace

def optimized_dna_to_dna_seq(dna_sequence, s_value,selected_host_organism):
    aa_sequence = dna_to_aa(dna_sequence)
    optimized_sequence,modify_seq_nospace = optimize_gene_sequence(aa_sequence, s_value,selected_host_organism)
    return optimized_sequence,modify_seq_nospace
