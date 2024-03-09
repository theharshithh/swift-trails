import os
from dotenv import dotenv_values
import streamlit as st
from PyPDF2 import PdfReader
from io import BytesIO
from openai import OpenAI

secrets = dotenv_values('.env')

openai_api_key = secrets["OPENAI_API_KEY"]

if not openai_api_key:
    print("OPENAI_API_KEY env variable is not set.")
else:
    pass

def read_pdf(uploaded_file):
    if uploaded_file is not None:
        pdf_content = uploaded_file.read()
        pdf_stream = BytesIO(pdf_content)  # Wrap the bytes content in a BytesIO object
        pdf_reader = PdfReader(pdf_stream)
        
        text_content = ""
        for page in pdf_reader.pages:
            text_content += page.extract_text()

        return text_content
    return None

def get_patient_condition(medical_report):
    client = OpenAI(api_key=openai_api_key)
    instructions = """I want u to act like a simple flagging bot. Lets keep it simple. I will share with you a medical report or a patient's condition, and I need u to give me the following thing: True if the patient is fit for the medical trials, false if the patient is not fit for the medical trials.
    Here is a reference of Vital Data: Patient ID,Age,Sex,Temperature,Heart Beat,Respiratory Rate,Systolic BP,Diastolic BP,Oxygen Saturation\n1,29,M,98.6,72,16,120,80,98\n2,34,F,96.9,75,14,115,75,99\n3,45,M,98.4,70,18,110,70,97\n4,26,F,97.8,78,16,118,76,98\n5,37,M,97.1,72,15,120,80,96\n6,50,F,98.6,74,16,122,82,97\n7,32,M,97.7,68,14,115,78,99\n8,28,F,98.5,76,18,117,75,98\n9,41,M,98.6,70,17,119,79,96\n10,28,M,97.06,72,16,120,80,98\n
        These are considered baseline data. If I give u a medical report/patient's health vitals. For it being true, the input values should be within the range of the baseline data. If the input values are not within the range, then it should be flagged as FALSE(if false, tell me why false- KEEP IT SHORT). If the input values are within the range, but not exactly the same, then it should be flagged as NEUTRAL."""

    response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {
                    "role": "system",
                    "content": instructions
                },
                {
                    "role": "user",
                    "content": ""
                },
                {
                    "role": "assistant",
                    "content": "Okay, I will give u True or False"
                },
                {
                    "role": "user",
                    "content": f"medical report - {medical_report}"
                },
            ],
            temperature=1,
            max_tokens=1000,
            top_p=1,
            frequency_penalty=0,
            presence_penalty=0
        )
    response_message = response.choices[0].message.content
    
    return response_message

def get_flag_value(patient_result):
    if "True" in patient_result or "eligible" in patient_result.lower():
        return "Eligible"
    else:
        return "Ineligible"

def add_patient(patient_name, age, medical_report, flag_result, patient_result):
    # Use Streamlit's session state to store the reviewed_patients list
    if "reviewed_patients" not in st.session_state:
        st.session_state.reviewed_patients = []

    st.session_state.reviewed_patients.append({"Patient Name": patient_name, "Patient Age": age, "Flagging Result": flag_result, "Explanation": patient_result})

    # Display the reviewed patient details in a table below the page
    st.write("Reviewed Patient Details:")
    st.table(st.session_state.reviewed_patients)