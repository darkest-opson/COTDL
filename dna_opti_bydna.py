from colorama import Fore, Style
import streamlit as st
from full_opti_dict import *

def is_dna_sequence(sequence):
        valid_nucleotides = {'A', 'T', 'G', 'C'}
        return all(nucleotide in valid_nucleotides for nucleotide in sequence)

def optimize_codons(dna_seq, selected_host_organism):
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    if selected_host_organism =='Cicer arientum'or selected_host_organism=='Lathyrus sativus' or selected_host_organism=='Vigna mungo':
        codon_bias=codon_bias_cvml
    if selected_host_organism == 'Cajanus cajan' or selected_host_organism=='Phaseolus vulgaris' or selected_host_organism=='Vigna radiata'or selected_host_organism == 'Glycine max' :
        codon_bias=codon_bias_cgpv
    if selected_host_organism == 'Pisum sativum' or selected_host_organism=='lens culinaris':
        codon_bias=codon_bias_pslc

    optimized_seq = ''
    optimized_seq_ori = ''

    # Divide the DNA sequence into codons and replace each codon with the preferred codon
    for i in range(0, len(dna_seq), 3):
        if not is_dna_sequence(dna_seq):
            print("Input sequence is NOT a valid DNA sequence")
            return None
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            aa = codon_table.get(codon, '')
            if aa:
                preferred_codon = codon_bias.get(aa, codon)  # Get preferred codon or keep original if not specified
                if preferred_codon != codon:
                    optimized_seq += preferred_codon
                    optimized_seq_ori += preferred_codon
                else:
                    optimized_seq += codon
                    optimized_seq_ori += codon
            else:
                optimized_seq += codon
                optimized_seq_ori += codon

    return optimized_seq, optimized_seq_ori

# # Example usage:
# dna_sequence = input("Please enter the DNA sequence: ")
# formatted_sequence = ' '.join([dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)])
# st.write("Original sequence:", formatted_sequence)

# # Example codon bias
# codon_bias = {'F': 'TTT', 'L': 'CTT', 'I': 'ATT', 'M': 'ATG', 'V': 'GTT', 'S': 'TCT', 'P': 'CCT', 'T': 'ACT', 'A': 'GCT', 'Y': 'TAT', 'H': 'CAT', 'Q': 'CAA', 'N': 'AAT', 'K': 'AAA', 'D': 'GAT', 'E': 'GAA', 'C': 'TGT', 'W': 'TGG', 'R': 'AGA', 'G': 'GGA'}
# optimized_sequence = optimize_codons(dna_sequence)

# st.write("Optimized sequence:")
# st.write(optimized_sequence)
