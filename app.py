import streamlit as st
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
# from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
import requests
import streamlit as st
from PyPDF2 import PdfReader
from io import BytesIO
from openai import OpenAI

#file imports
from data_management import get_classyfire_info
from tox_as_api import toxicity_organ_mapping, get_toxicity_data, options 
from patient_manage import read_pdf, add_patient, get_patient_condition, get_flag_value
from complete_tox_report import get_parameters

flag_result ='' 
reviewed_patients = []

# Home page dir; 
def home_page():
    st.title("Swift-Trials v1.4")
    st.markdown("<div class='page-container'>", unsafe_allow_html=True)
    st.write("#### Drug Data:")
    smiles_input_key = "smiles_input"  # Unique key for the input field
    smiles = st.text_input("Enter the SMILES string", key=smiles_input_key, placeholder="Eg. CC(=O)OC1=CC=CC=C1C(=O)O")
    if st.button("Search", key='search_btn_home'):
        inchi_key, superclass, class_, subclass, molecular_framework, pathway = get_classyfire_info(smiles)
        compound = Chem.MolFromSmiles(smiles)
        if compound is not None:
            # name = compound.GetProp("_Name") if "_Name" in compound.GetPropNames() else "N/A"
            formula = rdMolDescriptors.CalcMolFormula(compound)
            weight = rdMolDescriptors.CalcExactMolWt(compound)
                
            st.write("Molecular Formula:", formula)
            st.write("Molecular Weight:", str(round(weight, 2)))
            st.write(f"InChIKey: {inchi_key}")
            st.write(f"Class: {class_}")
            st.write(f"Subclass: {subclass}")
            st.write(f"Superclass: {superclass}")
            st.write(f"Molecular Framework: {molecular_framework}")
        else:
            st.write("Invalid compound generated from SMILES.")

def toxicity():
    st.title("Toxicity Report")
    st.markdown("<div class='page-container'>", unsafe_allow_html=True)
    smile_st = st.text_input("Search for Individual Toxicity Parameter", placeholder="SMILE String Eg. CCCN")
    toxicity_param = st.selectbox('Select a Parameter:', options)
    if st.button("Submit"):
        st.write(get_toxicity_data(smile_st, toxicity_param))
    st.markdown("</div>", unsafe_allow_html=True)
    complete_report_btn = st.button("Access Complete Report", key="center_aligned_button")
    if(complete_report_btn):
        parameters = get_parameters(smile_st)
        if isinstance(parameters, dict):
            st.write("Properties for the provided SMILES structure:")
            st.table(pd.DataFrame(parameters))
        else:
            st.write(parameters)


def patient_flagging_page():
    st.title("Patient Management")
    st.write("Flagging Patients for Clinical Trails")

    # Collect name and age as inputs
    name = st.text_input('Enter Patient Name:', placeholder='John')
    age = st.number_input('Enter Patient Age:', placeholder='31', max_value=150, step=1)
    uploaded_file = st.file_uploader("Choose a PDF file", type=["pdf"])

    if uploaded_file is not None:
        medical_report = read_pdf(uploaded_file)
        if medical_report is not None:
            st.write("PDF Content:")
            st.write(medical_report)
            if st.button('Submit'):
                result = get_patient_condition(medical_report)
                st.write(result)

                flag_val = get_flag_value(result)
                add_patient(name, age, medical_report, flag_val, result)
    else:
        st.write("Please upload a valid PDF file.")

# Navigation
pages = {
    "Home": home_page,
    "Toxicity Report": toxicity,
    "Patient Flagging": patient_flagging_page
    # , 
    # "Similarity Score": 
}

selection = st.sidebar.radio("Navigate:", list(pages.keys()))
pages[selection]()
