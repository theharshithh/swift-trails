import streamlit as st
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import requests

def get_compound_name_from_smiles(smiles):
    try:
        # Search PubChem with the provided SMILES
        compound = pcp.get_compounds(smiles, 'smiles', record_type='3d')[0]

        # Get the compound name
        compound_name = compound.iupac_name

        return compound_name

    except Exception as e:
        print(f"Error: {e}")
        return None

def get_classyfire_info(smiles):
    base_url = "https://structure.gnps2.org/classyfire"
    params = {"smiles": smiles}

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an exception for bad requests
        data = response.json()

        # Check if the response contains the expected keys
        if "inchikey" in data and "superclass" in data:
            inchl_key = data["inchikey"]
            superclass = data["superclass"].get("name", "")
            class_ = data.get("class", {}).get("name", "")
            
            # Check for the presence of 'subclass' key before accessing it
            subclass_data = data.get("subclass", {})
            subclass = subclass_data.get("name", "") if subclass_data else ""

            molecular_framework = data.get("molecular_framework", "")
            pathway = data.get("pathway", "")

            return inchl_key, superclass, class_, subclass, molecular_framework, pathway
        else:
            print("Error: Unexpected response format")
            return "", "", "", "", "", ""

    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")
        return "", "", "", "", "", ""
