import streamlit as st
from Bio import SeqIO
from collections import Counter
import io
from dna_opti_bydna import optimize_codons
from res_check import check_res_site_function
from Bio.Restriction import *
import sys
from io import StringIO
import pandas as pd
from io import BytesIO
from Bio.Restriction.PrintFormat import PrintFormat
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from prob_opti_dna import optimized_dna_to_dna_seq,optimize_gene_sequence
from full_opti_aa import optimize_dna_sequence
from gc_content import calculate_gc_content,calculate_melting_temp
from neatbio import sequtils as utils
from deoptimize import deoptimize_sequence_aa,deoptimized_dna_to_dna_seq
# data vis pkgs
import matplotlib.pyplot as plt
import matplotlib
import json
matplotlib.use("Agg")


def main():
    flag = False
    st.title("CODTL: Codon Optimization and Deoptimization Tool for Legume")

    menu = ["Home", "Optimization", "Deoptimization", "Tutorial"]
    choice = st.sidebar.selectbox("Select Activity", menu)

    if choice == "Home":
        st.markdown("""
        Welcome to CODTL, the Legume DNA Optimization Tool. This application enables users to optimize DNA sequences for high protein expression and transformation experiments in legumes. Upload DNA or amino acid sequences and apply codon bias optimization techniques for enhanced gene expression in legumes. Ideal for researchers, biologists, and genetic engineers, CODTL offers the essential tools for designing legume-specific DNA sequences.
        """)
        st.image("rf.jpg", use_column_width=True)


    elif choice == "Optimization":
  
        st.subheader("Get optimize DNA sequence.")

        if st.checkbox("By DNA sequence"):
            seq_file=st.file_uploader("Upload FASTA File",type=["fasta","fa","txt"])

            if seq_file is not None:
                # Convert the uploaded file object to a string buffer
                # print("Before reading the sequence file")  # Print statement added
                seq_buffer = io.StringIO(seq_file.getvalue().decode("utf-8"))
                # /print("After creating the string buffer")  # Print statement added
                dna_record = SeqIO.read(seq_buffer, "fasta")

                input_dna_seq=dna_record.seq
                modify_seq_ori = None



                details=st.radio("Details",("Check Description of Input file","Check Input Sequence"), index=None)
                if details=="Check Description of Input file":
                    st.write(dna_record.description)
                elif details=="Check Input Sequence":
                    st.write(dna_record.seq)

                host_organisms = ["Phaseolus vulgaris", "Lathyrus sativus", "Vigna mungo", "Vigna radiata","Cajanus cajan", "Glycine max", "Pisum sativum", "lens culinaris", "Cicer arientum"]
                selected_host_organism = st.selectbox("Select Host Organism", host_organisms)

                Get_output=st.radio("Get modify DNA Sequence",("Full Optimization","Probabilistic optimization"),index=None)


                if Get_output=="Full Optimization":
                    try:
                        flag = True
                        if len(input_dna_seq) % 3 != 0:
                            st.write("Error: DNA sequence length must be a multiple of 3")
                        else:
                            st.write("This might take few minutes.Please be patient.")
                            modify_seq,modify_seq_ori=optimize_codons(input_dna_seq)
                            st.markdown(f"Optimized sequence: {modify_seq}", unsafe_allow_html=True)

                            melting_temp = calculate_melting_temp(modify_seq)
                            gc_content = calculate_gc_content(modify_seq)
            
                            # Create a DataFrame for the Excel file
                            data = {
                                "Modified Sequence": [modify_seq],
                                "Length": [len(modify_seq)],
                                "Melting Temperature": [melting_temp],
                                "GC Content": [gc_content]
                            }
                            df = pd.DataFrame(data)
            
                            # Convert DataFrame to Excel
                            output = BytesIO()
                            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                                df.to_excel(writer, index=False, sheet_name='Sheet1')
                            processed_data = output.getvalue()

                            # Provide download link
                            st.download_button(
                                label="Download data as Excel",
                                data=processed_data,
                                file_name='optimized_dna_sequence.xlsx',
                                mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                            )
                    except:
                        flag = False
                        pass
                elif Get_output=="Probabilistic optimization":
                    try:
                        flag = True
                        st.write("In the probalistics approch value of Optimality factor(S) is taken\n 1.S=0 signifies full optimization\n 2.S=64 signifies anti-optimization\n 3.S=1 signifies codons in the genrated sequence randomly selected.")
                        # S=st.number_input("Enter The value of 'S':")
                        optimality_factor = st.number_input("Enter Optimality Factor (0 to 64)", min_value=0, max_value=64, step=1)
                        confirm = st.button("Confirm selection")
                        
                        modify_seq,modify_seq_ori=optimized_dna_to_dna_seq(input_dna_seq,optimality_factor,selected_host_organism)
                        st.markdown(f"Optimized sequence: {modify_seq}", unsafe_allow_html=True)

                        melting_temp = calculate_melting_temp(modify_seq)
                        gc_content = calculate_gc_content(modify_seq)
            
                            # Create a DataFrame for the Excel file
                        data = {
                                "Modified Sequence": [modify_seq],
                                "Length": [len(modify_seq)],
                                "Melting Temperature": [melting_temp],
                                "GC Content": [gc_content]
                            }
                        df = pd.DataFrame(data)
            
                            # Convert DataFrame to Excel
                        output = BytesIO()
                        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                            df.to_excel(writer, index=False, sheet_name='Sheet1')
                        processed_data = output.getvalue()

                            # Provide download link
                        st.download_button(
                                label="Download data as Excel",
                                data=processed_data,
                                file_name='optimized_dna_sequence.xlsx',
                                mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                            )
                    except:
                        st.write("Input sequence is NOT a valid DNA sequence")
                        flag = False
                        pass
                if flag == True:
                    check_res_site=st.radio("Check for the Restriction site in the DNA sequences",("Check RE in input sequence","Check RE in modified sequence"),index=None)
                    if check_res_site=="Check RE in input sequence":
                        ana=check_res_site_function(input_dna_seq)
                        st.write(str(ana.full()))
                        kj=ana.full()



                    elif check_res_site=="Check RE in modified sequence":
                        ana_1=check_res_site_function(modify_seq_ori)
                        st.write(str(ana_1.full()))

                    st.subheader("Nucleotide Frequency")
                    dna_freq = Counter(modify_seq_ori)
                    st.write(dna_freq)

                    if st.button("Plot Freq"):
                        fig, ax = plt.subplots(figsize=(5, 3))  # Adjust the figsize as per your preference
                        barlist = plt.bar(dna_freq.keys(), dna_freq.values())
                        st.pyplot(fig)

                    st.subheader("DNA Composition")
                    gc_score=utils.gc_content(str(modify_seq_ori))
                    at_score=utils.at_content(str(modify_seq_ori))
                    st.json({"GC Content":gc_score,"AT Content":at_score})

                    try:
                        # protein synthesis
                        st.subheader("Protein Sythesis")
                        modseq_1=Seq(modify_seq_ori)
                        pl=modseq_1.translate()
                        aa_freq = Counter(str(pl))

                        if st.checkbox("Transcribe Sequence of modified sequence"):
                            st.write(modseq_1.transcribe())
                        elif st.checkbox("Translated Sequence of modified sequence"):
                            st.write(modseq_1.translate())
                        elif st.checkbox("Complemented Sequence of modified sequence"):
                            st.write(modseq_1.complement())
                        elif st.checkbox("AA Frequency"):
                            st.write(aa_freq)
                        elif st.checkbox("Plot AA Frequency"):
                            fig, ax = plt.subplots(figsize=(5, 3))
                            barlist=plt.bar(aa_freq.keys(),aa_freq.values())
                            st.pyplot(fig)
                    except:
                        st.write("Choose an option first between full optimisation and probablistic optimisation")


                

        elif st.checkbox("By protein sequence"):
            flag = False
            seq_file=st.file_uploader("Upload FASTA File",type=["fasta","fa","txt"])

            if seq_file is not None:
                    # Convert the uploaded file object to a string buffer
                    # print("Before reading the sequence file")  # Print statement added
                seq_buffer = io.StringIO(seq_file.getvalue().decode("utf-8"))
                    # /print("After creating the string buffer")  # Print statement added
                aa_record = SeqIO.read(seq_buffer, "fasta")

                input_aa_seq=aa_record.seq
                modify_seq_ori_aa = None

                details=st.radio("Details",("Check Description of Input file","Check Input Sequence"), index=None)
                if details=="Check Description of Input file":
                    st.write(aa_record.description)
                elif details=="Check Input Sequence":
                    st.write(aa_record.seq)

                host_organisms = ["Phaseolus vulgaris", "Lathyrus sativus", "Vigna mungo", "Vigna radiata","Cajanus cajan", "Glycine max", "Pisum sativum", "lens culinaris", "Cicer arientum"]
                selected_host_organism = st.selectbox("Select Host Organism", host_organisms)

                Get_output=st.radio("Get Modified DNA Sequence",("Full Optimization","Probabilistic optimization"),index=None)
                if Get_output=="Full Optimization":
                    flag = True
                    st.write("This might take few minutes.Please be patient.")
                    modify_seq_aa,modify_seq_ori_aa=optimize_dna_sequence(input_aa_seq)
                    st.markdown(f"Optimized sequence: {modify_seq_aa}", unsafe_allow_html=True)
                    melting_temp = calculate_melting_temp(modify_seq_aa)
                    gc_content = calculate_gc_content(modify_seq_aa)
            
                        # Create a DataFrame for the Excel file
                    data = {
                            "Modified Sequence": [modify_seq_aa],
                            "Length": [len(modify_seq_aa)],
                            "Melting Temperature": [melting_temp],
                            "GC Content": [gc_content]
                        }
                    df = pd.DataFrame(data)
            
                        # Convert DataFrame to Excel
                    output = BytesIO()
                    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                        df.to_excel(writer, index=False, sheet_name='Sheet1')
                    processed_data = output.getvalue()

                        # Provide download link
                    st.download_button(
                            label="Download data as Excel",
                            data=processed_data,
                            file_name='optimized_dna_sequence.xlsx',
                            mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                        )

                elif Get_output=="Probabilistic optimization":
                    flag = True
                    st.write("In the probalistics approch value of Optimality factor(S) is taken\n 1.S=0 signifies full optimization\n 2.S=64 signifies anti-optimization\n 3.S=1 signifies codons in the genrated sequence randomly selected.")
                        # S=st.number_input("Enter The value of 'S':")
                    optimality_factor = st.number_input("Enter Optimality Factor (0 to 64)", min_value=0, max_value=64, step=1)
                    confirm = st.button("Confirm selection")
                    # if confirm:
                    modify_seq_aa,modify_seq_ori_aa=optimize_gene_sequence(input_aa_seq,optimality_factor,selected_host_organism)
                    st.markdown(f"Optimized sequence: {modify_seq_aa}", unsafe_allow_html=True)

                    melting_temp = calculate_melting_temp(modify_seq_aa)
                    gc_content = calculate_gc_content(modify_seq_aa)
            
                            # Create a DataFrame for the Excel file
                    data = {
                                "Modified Sequence": [modify_seq_aa],
                                "Length": [len(modify_seq_aa)],
                                "Melting Temperature": [melting_temp],
                                "GC Content": [gc_content]
                            }
                    df = pd.DataFrame(data)
            
                            # Convert DataFrame to Excel
                    output = BytesIO()
                    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                        df.to_excel(writer, index=False, sheet_name='Sheet1')
                    processed_data = output.getvalue()

                            # Provide download link
                    st.download_button(
                                label="Download data as Excel",
                                data=processed_data,
                                file_name='optimized_dna_sequence.xlsx',
                                mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                            )

                if flag == True:                        
                    try:
                        st.subheader("Check Restriction Site in optimized DNA Seq")
                        if st.checkbox("Check Restriction Site in optimized DNA Seq"):
                            mod=Seq(modify_seq_ori_aa)
                            ana_1=check_res_site_function(mod)
                            st.write(str(ana_1.full()))
                    except Exception as e:
                        st.write(f"Choose Full/Probabilistic optimisation or Sequence may contains stop codon.")
                        st.write("[ '*'  in nucleotide frequency represents stop codon.]")


                    st.subheader("Nucleotide Frequency")
                    dna_freq = Counter(modify_seq_ori_aa)
                    st.write(dna_freq)

                    if st.button("Plot Freq"):
                        fig, ax = plt.subplots(figsize=(5, 3))  # Adjust the figsize as per your preference
                        barlist = plt.bar(dna_freq.keys(), dna_freq.values())
                        st.pyplot(fig)

                    st.subheader("DNA Composition")
                    gc_score=utils.gc_content(str(modify_seq_ori_aa))
                    at_score=utils.at_content(str(modify_seq_ori_aa))
                    st.json({"GC Content":gc_score,"AT Content":at_score})

                    # protein synthesis
                    try:
                        st.subheader("Protein Sythesis")
                        modseq=Seq(modify_seq_ori_aa)
                        # print(modseq)
                        pl=modseq.translate()
                        aa_freq = Counter(str(pl))

                        if st.checkbox("Transcribe Sequence of modified sequence"):
                            st.write(modseq.transcribe())
                        elif st.checkbox("Translated Sequence of modified sequence"):
                            st.write(modseq.translate())
                        elif st.checkbox("Complemented Sequence of modified sequence"):
                            st.write(modseq.complement())
                        elif st.checkbox("AA Frequency"):
                            st.write(aa_freq)
                        elif st.checkbox("Plot AA Frequency"):
                            fig, ax = plt.subplots(figsize=(5, 3))
                            barlist=plt.bar(aa_freq.keys(),aa_freq.values())
                            st.pyplot(fig)
                    except Exception as e :
                            st.write(f"Choose Full/Probabilistic optimisation or Sequence may contains stop codon.")
                            st.write("[ '*'  in nucleotide frequency represents stop codon.]")

                

                




    elif choice == "Deoptimization":
        try:
            st.markdown("""
            Choose specific amino acid for Deoptimization.
            """)
            details_1 = st.radio("Input File is a:", ("By DNA sequence", "By amino acid sequence"), index=None)
            if details_1 == "By DNA sequence":
                st.subheader("Get deoptimize DNA sequence by input DNA sequence.")
                seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa", "txt"])

                if seq_file is not None:
                    # Convert the uploaded file object to a string buffer
                    # print("Before reading the sequence file")  # Print statement added
                    seq_buffer = io.StringIO(seq_file.getvalue().decode("utf-8"))
                    # /print("After creating the string buffer")  # Print statement added
                    dna_record = SeqIO.read(seq_buffer, "fasta")

                    input_dna_seq = dna_record.seq
                    demodify_ori = None

                    details = st.radio("Details", ("Check Description of Input file", "Check Input Sequence"), index=None)
                    if details == "Check Description of Input file":
                        st.write(dna_record.description)
                    elif details == "Check Input Sequence":
                        st.write(dna_record.seq)


                if st.checkbox("Deoptimized Sequence"):

                    host_organisms = ["Phaseolus vulgaris", "Lathyrus sativus", "Vigna mungo", "Vigna radiata","Cajanus cajan", "Glycine max", "Pisum sativum", "lens culinaris", "Cicer arientum"]
                    selected_host_organism = st.selectbox("Select Host Organism", host_organisms)

                    amino_acid_names = {
                        'M': 'Methionine',
                        'F': 'Phenylalanine',
                        'L': 'Leucine',
                        'S': 'Serine',
                        'Y': 'Tyrosine',
                        'C': 'Cysteine',
                        'W': 'Tryptophan',
                        'P': 'Proline',
                        'Q': 'Glutamine',
                        'H': 'Histidine',
                        'R': 'Arginine',
                        'I': 'Isoleucine',
                        'T': 'Threonine',
                        'N': 'Asparagine',
                        'K': 'Lysine',
                        'D': 'Aspartic acid',
                        'E': 'Glutamic acid',
                        'V': 'Valine',
                        'A': 'Alanine',
                        'G': 'Glycine'
                    }
                    amino_acids = ['M', 'F', 'L', 'S', 'Y', 'C', 'W', 'P', 'Q', 'H', 'R', 'I', 'T', 'N', 'K', 'D', 'E', 'V', 'A', 'G']
                    selected_amino_acids = st.multiselect('Choose amino acids to deoptimize:',
                                                        amino_acids,
                                                        format_func=lambda x: f"{x}:{amino_acid_names[x]}")

                    if st.checkbox("Get deoptimize DNA seq for selected amino acid"):
                        de_opt,demodify_ori=deoptimized_dna_to_dna_seq(input_dna_seq,selected_amino_acids,selected_host_organism)
                        st.markdown(f"deoptimized sequence: {de_opt}", unsafe_allow_html=True)

                        
                        melting_temp = calculate_melting_temp(demodify_ori)
                        gc_content = calculate_gc_content(demodify_ori)
            
                        # Create a DataFrame for the Excel file
                        data = {
                            "Modified Sequence": [demodify_ori],
                            "Length": [len(demodify_ori)],
                            "Melting Temperature": [melting_temp],
                            "GC Content": [gc_content]
                        }
                        df = pd.DataFrame(data)
            
                        # Convert DataFrame to Excel
                        output = BytesIO()
                        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                            df.to_excel(writer, index=False, sheet_name='Sheet1')
                        processed_data = output.getvalue()

                        # Provide download link
                        st.download_button(
                            label="Download data as Excel",
                            data=processed_data,
                            file_name='optimized_dna_sequence.xlsx',
                            mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                        )


            elif details_1 == "By amino acid sequence":
                st.subheader("Get deoptimize DNA sequence by input Amino acid sequence.")
                seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa", "txt"])

                if seq_file is not None:
                # Convert the uploaded file object to a string buffer
                # print("Before reading the sequence file")  # Print statement added
                    seq_buffer = io.StringIO(seq_file.getvalue().decode("utf-8"))
                    # /print("After creating the string buffer")  # Print statement added
                    aa_record = SeqIO.read(seq_buffer, "fasta")

                    input_aa_seq=aa_record.seq
                    modify_seq_ori_aa = None

                    details=st.radio("Details",("Check Description of Input file","Check Input Sequence"), index=None)
                    if details=="Check Description of Input file":
                        st.write(aa_record.description)
                    elif details=="Check Input Sequence":
                        st.write(aa_record.seq)

                if st.checkbox("Deoptimized Sequence"):
                    amino_acid_names = {
                        'M': 'Methionine',
                        'F': 'Phenylalanine',
                        'L': 'Leucine',
                        'S': 'Serine',
                        'Y': 'Tyrosine',
                        'C': 'Cysteine',
                        'W': 'Tryptophan',
                        'P': 'Proline',
                        'Q': 'Glutamine',
                        'H': 'Histidine',
                        'R': 'Arginine',
                        'I': 'Isoleucine',
                        'T': 'Threonine',
                        'N': 'Asparagine',
                        'K': 'Lysine',
                        'D': 'Aspartic acid',
                        'E': 'Glutamic acid',
                        'V': 'Valine',
                        'A': 'Alanine',
                        'G': 'Glycine'
                    }
                    
                    host_organisms = ["Phaseolus vulgaris", "Lathyrus sativus", "Vigna mungo", "Vigna radiata","Cajanus cajan", "Glycine max", "Pisum sativum", "lens culinaris", "Cicer arientum"]
                    selected_host_organism = st.selectbox("Select Host Organism", host_organisms)
                    
                    amino_acids = ['M', 'F', 'L', 'S', 'Y', 'C', 'W', 'P', 'Q', 'H', 'R', 'I', 'T', 'N', 'K', 'D', 'E', 'V', 'A', 'G']
                    selected_amino_acids = st.multiselect('Choose amino acids to deoptimize:',
                                                        amino_acids,
                                                        format_func=lambda x: f"{x}:{amino_acid_names[x]}")
                    

                    if st.checkbox("Get deoptimize DNA seq for selected amino acid"):
                        de_opt,demodify_ori=deoptimize_sequence_aa(input_aa_seq,selected_amino_acids,selected_host_organism)
                        st.markdown(f"deoptimized sequence: {de_opt}", unsafe_allow_html=True)

                        melting_temp = calculate_melting_temp(demodify_ori)
                        gc_content = calculate_gc_content(demodify_ori)
            
                        # Create a DataFrame for the Excel file
                        data = {
                            "Modified Sequence": [demodify_ori],
                            "Length": [len(demodify_ori)],
                            "Melting Temperature": [melting_temp],
                            "GC Content": [gc_content]
                        }
                        df = pd.DataFrame(data)
            
                        # Convert DataFrame to Excel
                        output = BytesIO()
                        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                            df.to_excel(writer, index=False, sheet_name='Sheet1')
                        processed_data = output.getvalue()

                        # Provide download link
                        st.download_button(
                            label="Download data as Excel",
                            data=processed_data,
                            file_name='optimized_dna_sequence.xlsx',
                            mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                        )
        except:
            st.write("Please provide a valid input sequence.")


        

        

    
    elif choice == "Tutorial":
        st.subheader("Optimization")

        st.markdown("""
        ### Step 1: Choose Sequence Type
        In the optimization section, choose "DNA sequence" if your input file is a DNA sequence, or choose "Amino acid sequence" if your input file contains protein amino acids.

        ### Step 2: Upload Sequence File
        Upload your sequence file. Ensure your file is correctly formatted to avoid any errors during the optimization process.

        ### Step 3: Select Host Organism
        Select your desired expression host organism from the dropdown menu. This step is crucial as different organisms have different codon preferences.

        ### Step 4: Optimize Your Sequence
        Get your optimized DNA by clicking either "Full Optimization" or "Probabilistic Optimization." For probabilistic optimization, you need to choose the Optimality factor (S):
        - **S = 0**: Signifies full optimization.
        - **S = 64**: Signifies anti-optimization.
        - **S = 1**: Codons in the generated sequence are randomly selected.


        ### Step 5: Download Optimized Sequence
        You can download the optimized sequence in Excel format for further analysis and use.

        ### Additional Features
        - **Check Restriction Sites:** Analyze the restriction sites in both the input and optimized sequences.
        - **GC Content:** Calculate the GC content of your sequences.
        - **Transcription and Translation:** Get the transcribed and translated sequence of the optimized DNA.
        """)

        st.markdown("""
        ### Deoptimization Tutorial

        #### Step 1: Choose Specific Amino Acid for Deoptimization
        Select the amino acid(s) you wish to deoptimize.

        #### Step 2: Choose Sequence Type
        Choose "DNA sequence" if your input file is a DNA sequence, or "Amino acid sequence" if your input file contains protein amino acids.

        #### Step 3: Upload Sequence File
        Upload your sequence file. Ensure your file is correctly formatted to avoid any errors during the deoptimization process.

        #### Step 4: Select Host Organism
        Select your desired expression host organism from the dropdown menu. This step is crucial as different organisms have different codon preferences.

        #### Step 5: Choose Amino Acids to Deoptimize
        Select the amino acids you wish to deoptimize from the list provided.

        #### Step 6: Get Deoptimized Sequence
        Check the checkbox "Get deoptimized DNA sequence for selected amino acid" to generate the deoptimized sequence.

        #### Step 7: Download Deoptimized Sequence
        You can download the deoptimized sequence in Excel format for further analysis and use.
        """)



if __name__ == "__main__":
    main()
