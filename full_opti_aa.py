# codon_bias = {
#     'A': 'GCT',  # Alanine
#     'C': 'TGT',  # Cysteine
#     'D': 'GAT',  # Aspartic acid
#     'E': 'GAA',  # Glutamic acid
#     'F': 'TTT',  # Phenylalanine
#     'G': 'GGA',  # Glycine
#     'H': 'CAT',  # Histidine
#     'I': 'ATT',  # Isoleucine
#     'K': 'AAA',  # Lysine
#     'L': 'CTT',  # Leucine
#     'M': 'ATG',
#     'N': 'AAT',  # Asparagine
#     'P': 'CCT',  # Proline
#     'Q': 'CAA',  # Glutamine
#     'R': 'AGA',  # Arginine
#     'S': 'TCT',  # Serine
#     'T': 'ACT',  # Threonine
#     'V': 'GTT',  # Valine
#     'Y': 'TAT'   # Tyrosine
    
# }
def optimize_dna_sequence(aa_sequence):

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
    'M': 'ATG',
    'N': 'AAT',  # Asparagine
    'P': 'CCT',  # Proline
    'Q': 'CAA',  # Glutamine
    'R': 'AGA',  # Arginine
    'S': 'TCT',  # Serine
    'T': 'ACT',  # Threonine
    'V': 'GTT',  # Valine
    'Y': 'TAT'   # Tyrosine
    
    }   
    # Initialize an empty string to store the optimized DNA sequence
    dna_sequence = ''
    opti_seq=''
    
    # Iterate over each amino acid in the amino acid sequence
    for aa in aa_sequence:
        # Lookup the preferred codon for the current amino acid in the codon bias table
        preferred_codon = codon_bias.get(aa, '')
        
        # If a preferred codon is found, append it to the DNA sequence
        if preferred_codon:
            dna_sequence += preferred_codon
            opti_seq+=preferred_codon
        else:
            # If no preferred codon is found (e.g., for stop codons), append an empty string
            dna_sequence += '*'
            opti_seq+='*'
    
    return dna_sequence,opti_seq
