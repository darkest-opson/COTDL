from codon_dict import *
def is_dna_sequence(sequence):
        valid_nucleotides = {'A', 'T', 'G', 'C'}
        return all(nucleotide in valid_nucleotides for nucleotide in sequence)

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




def deoptimize_sequence_aa(seq, selected_aa,selected_host_organism):
    codon_bias = {
        'A': 'GCT',  # Alanine
        'C': 'TGT',  # Cysteine
        'D': 'GAT',  # Aspartic acid
        'E': 'GAA',  # Glutamic acid
        'F': 'TTT',  # Phenylalanine
        'G': 'GGA',  # Glycine
        'H': 'CAT',  # Histidine
        'I': 'ATT',  # Isoleucine
        'K': 'AAA',  # Lysine
        'L': 'CTT',  # Leucine
        'M': 'ATG',   # Methionine
        'N': 'AAT',  # Asparagine
        'P': 'CCT',  # Proline
        'Q': 'CAA',  # Glutamine
        'R': 'AGA',  # Arginine
        'S': 'TCT',  # Serine
        'T': 'ACT',  # Threonine
        'V': 'GTT',  # Valine
        'Y': 'TAT'   # Tyrosine
    }
    
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

    deoptimized_dna_sequence = ""
    deopti_seq = ""
    
    for aa in seq:
        if aa in selected_aa:
            codon_frequencies = codon_bias_dict[aa]
            least_frequent_codon = min(codon_frequencies, key=codon_frequencies.get)
            deoptimized_dna_sequence += least_frequent_codon + " "
            deopti_seq += least_frequent_codon
        else:
            # Lookup the preferred codon for the current amino acid in the codon bias table
            preferred_codon = codon_bias.get(aa, '')
            
            # If a preferred codon is found, append it to the DNA sequence
            if preferred_codon:
                deoptimized_dna_sequence += preferred_codon + " "
                deopti_seq += preferred_codon
            else:
                # If no preferred codon is found (e.g., for stop codons), append an empty string
                deoptimized_dna_sequence += '*' + " "
                deopti_seq += '*'
    
    return deoptimized_dna_sequence, deopti_seq

# def deoptimized_dna_to_dna_seq(dna_sequence,selected_aa):
#     aa_sequence = dna_to_aa(dna_sequence)
#     deoptimized_sequence,demodify_seq_nospace = deoptimize_sequence_aa(aa_sequence,selected_aa)
#     return deoptimized_sequence,demodify_seq_nospace
# ******************************************************************************************************************

def translate_codon(codon,selected_host_organism):
    """Translate a codon to its corresponding amino acid using the codon bias dictionary."""
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
        
    for aa, codons in codon_bias_dict.items():
        if codon in codons:
            return aa
    return None  # Return None if codon is not found in the dictionary


def is_dna_sequence(seq):
    """Check if a sequence is a valid DNA sequence."""
    valid_bases = set('ATGC')
    return all(base in valid_bases for base in seq.upper())



def deoptimized_dna_to_dna_seq(dna_seq, selected_aa, selected_host_organism):
    if not is_dna_sequence(dna_seq):
        st.write("Input sequence is NOT a valid DNA sequence")
        return None
    
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
        
    deoptimized_dna_sequence = ""
    deopti_seq = ""
    
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            # Translate codon to amino acid
            aa = translate_codon(codon,selected_host_organism)
            if aa in selected_aa:
                # Lookup the preferred codon for the current amino acid in the codon bias table
                codon_frequencies = codon_bias_dict[aa]
                least_frequent_codon = min(codon_frequencies, key=codon_frequencies.get)
                deoptimized_dna_sequence += f"<span style='color:red'>{least_frequent_codon}</span>"
                deopti_seq += least_frequent_codon
            else:
                # Keep the sequence unchanged if the amino acid is not selected for deoptimization
                deoptimized_dna_sequence += codon
                deopti_seq += codon
    return deoptimized_dna_sequence, deopti_seq



