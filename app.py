import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd
import io

def analyse_seq(data_1, igg1_constant, light_constant):
    seq_list = []
    seq_list_1 = []
    
    heavy_strings = io.StringIO()
    light_strings = io.StringIO()
    
    for i, row in data_1.iterrows():
        if pd.notna(row[2]):
            seq_list.append(row[1])
            heavy_strings.write(row[0]+'_'+row[1] + '\n')
    
    for i, row in data_1.iterrows():
        if pd.notna(row[1]):
            seq_list_1.append(row[2])
            light_strings.write(row[2] + '\n')
    
    heavy_strings.seek(0)
    light_strings.seek(0)
    
    result_list = []
    
    for record in SeqIO.parse(io.StringIO(igg1_constant), "fasta"):
        sequence_data = str(record.seq)
        header = record.description
        
        for record_light in SeqIO.parse(io.StringIO(light_constant), "fasta"):
            sequence_data_light = str(record_light.seq)
            
            strings_heavy = heavy_strings.readlines()
            strings_light = light_strings.readlines()
            
            for string_final_heavy, string_final_light in zip(strings_heavy, strings_light):    
                if string_final_heavy.startswith('>'):
                    seq_h_heavy = string_final_heavy.strip().split('_')
                    seq_hh_heavy = seq_h_heavy[1]
                    ses_index_heavy = seq_h_heavy[0]
                    seq_h_light = string_final_light.strip()
                    total_string = (seq_hh_heavy + sequence_data + seq_h_light + sequence_data_light)
                    
                    protein_analysis = ProteinAnalysis(total_string)
                    isolectric_point = protein_analysis.isoelectric_point()
                    mw = protein_analysis.molecular_weight()
                    extinction_coeff = protein_analysis.molar_extinction_coefficient()
                    abs_ec = round(extinction_coeff[1] / mw, 3)
                    result_list.append([ses_index_heavy, seq_hh_heavy, seq_h_light, isolectric_point, abs_ec, total_string])
    
    result_df = pd.DataFrame(result_list, columns=['Index', 'Heavy Chain', 'Light Chain', 'Pi', 'Abs EC', 'Total_string'])
    return result_df

st.title('Protein Sequence Analysis App')

st.write("Please upload the required files:")

data_1_file = st.file_uploader("Upload CSV file with sequences", type="csv")
igg1_constant_file = st.file_uploader("Upload IgG1 constant file", type="txt")
light_constant_file = st.file_uploader("Upload light constant file", type="txt")

if data_1_file and igg1_constant_file and light_constant_file:
    data_1 = pd.read_csv(data_1_file)
    igg1_constant = igg1_constant_file.getvalue().decode("utf-8")
    light_constant = light_constant_file.getvalue().decode("utf-8")
    
    if st.button('Analyze Sequences'):
        result_df = analyse_seq(data_1, igg1_constant, light_constant)
        st.write("Analysis complete. Here's a preview of the results:")
        st.write(result_df.head())
        
        csv = result_df.to_csv(index=False)
        st.download_button(
            label="Download results as CSV",
            data=csv,
            file_name="protein_analysis_results.csv",
            mime="text/csv",
        )
else:
    st.write("Please upload all required files to proceed with the analysis.")
