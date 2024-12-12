import requests
import pandas as pd
import streamlit as st
from pubchempy import get_compounds
from rdkit import Chem
import urllib.parse

# Helper: Validate and Canonicalize SMILES
def validate_and_canonicalize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            st.error("Invalid SMILES string. Could not be parsed.")
            return None
    except Exception as e:
        st.error(f"Error during SMILES validation: {e}")
        return None

# Fetch ChEMBL Data
def fetch_chembl_data(smiles):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?filter=molecule_structures__canonical_smiles__exact={urllib.parse.quote(smiles)}"
    headers = {"Accept": "application/json"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        try:
            data = response.json()
            molecules = [
                {
                    "ChEMBL ID": molecule.get("molecule_chembl_id", "N/A"),
                    "Name": molecule.get("pref_name", "N/A"),
                    "Max Phase": molecule.get("max_phase", "N/A"),
                }
                for molecule in data.get("molecules", [])
            ]
            return molecules
        except ValueError as e:
            st.error(f"JSON decode error: {e}")
            st.write("Raw Response Content:", response.text)
            return []
    else:
        st.error(f"ChEMBL request failed. HTTP Status: {response.status_code}")
        return []

# Fetch PubChem Data
def fetch_pubchem_data(smiles):
    compounds = get_compounds(smiles, 'smiles')
    data = []
    for compound in compounds:
        data.append({
            "CID": compound.cid,
            "Molecular Weight": compound.molecular_weight,
            "Canonical SMILES": compound.canonical_smiles,
            "IUPAC Name": compound.iupac_name,
        })
    return data

# Fetch UniProt Data
def fetch_uniprot_data(protein_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={protein_name}&format=json"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            return response.json().get("results", [])
        except ValueError as e:
            st.error(f"Error decoding UniProt JSON: {e}")
            return []
    else:
        st.error(f"UniProt request failed. HTTP Status: {response.status_code}")
        return []

# Save Results to Excel
def save_to_excel(data, filename="results.xlsx"):
    output = BytesIO()
    try:
        df = pd.DataFrame(data)
        df.to_excel(output, index=False, engine='openpyxl')
        output.seek(0)
        return output
    except Exception as e:
        st.error(f"Error saving to Excel: {e}")
        return None

# Main Streamlit App
def main():
    st.title("Comprehensive Compound Analysis Tool")
    st.write("Analyze compounds using ChEMBL, PubChem, and UniProt APIs.")

    # Input: SMILES string
    smiles = st.text_input("Enter the SMILES string of your compound:")

    if st.button("Analyze Compound"):
        if smiles:
            # Validate SMILES
            canonical_smiles = validate_and_canonicalize_smiles(smiles)
            if not canonical_smiles:
                st.error("Invalid SMILES. Analysis cannot proceed.")
                return

            results = {}

            # Fetch ChEMBL data
            st.spinner("Fetching ChEMBL data...")
            chembl_data = fetch_chembl_data(canonical_smiles)
            results['ChEMBL'] = chembl_data if chembl_data else "No data found."
            st.subheader("ChEMBL Data")
            st.dataframe(pd.DataFrame(chembl_data))

            # Fetch PubChem data
            st.spinner("Fetching PubChem data...")
            pubchem_data = fetch_pubchem_data(canonical_smiles)
            results['PubChem'] = pubchem_data if pubchem_data else "No data found."
            st.subheader("PubChem Data")
            st.dataframe(pd.DataFrame(pubchem_data))

            # Fetch UniProt data (example target: P53)
            st.spinner("Fetching UniProt data for P53...")
            uniprot_data = fetch_uniprot_data("P53")
            results['UniProt'] = uniprot_data if uniprot_data else "No data found."
            st.subheader("UniProt Data")
            st.dataframe(pd.DataFrame(uniprot_data))

            # Save results to Excel
            excel_file = save_to_excel(results)
            if excel_file:
                st.download_button(
                    label="Download Results as Excel",
                    data=excel_file,
                    file_name="compound_analysis_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )

if __name__ == "__main__":
    main()


