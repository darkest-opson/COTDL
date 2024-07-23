from full_opti_dict import *
def optimize_dna_sequence(aa_sequence,selected_host_organism):

    if selected_host_organism =='Cicer arientum'or selected_host_organism=='Lathyrus sativus' or selected_host_organism=='Vigna mungo':
        codon_bias=codon_bias_cvml
    if selected_host_organism == 'Cajanus cajan' or selected_host_organism=='Phaseolus vulgaris' or selected_host_organism=='Vigna radiata'or selected_host_organism == 'Glycine max' :
        codon_bias=codon_bias_cgpv
    if selected_host_organism == 'Pisum sativum' or selected_host_organism=='lens culinaris':
        codon_bias=codon_bias_pslc


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
