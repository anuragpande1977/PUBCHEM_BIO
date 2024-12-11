import requests
import pandas as pd
import streamlit as st
import urllib.parse
from io import BytesIO

# Helper: Validate SMILES
def validate_smiles(smiles):
    """
    Basic validation for SMILES string. Ensure it's non-empty and well-formed.
    """
    return bool(smiles.strip())

# Helper: Convert SMILES to PubChem CID
def get_cid_from_smiles(smiles):
    """
    Convert a SMILES string to a PubChem CID.
    """
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"
    params = {'smiles': smiles}
    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID'][0]
        else:
            st.warning("No CID found for the given SMILES string.")
            return None
    else:
        st.error(f"Error fetching CID from PubChem. Status code: {response.status_code}")
        return None

# Helper: Fetch PubChem Compound Summary
def fetch_pubchem_summary(cid):
    """
    Fetch summary information for a compound from PubChem.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,MolecularFormula,IUPACName/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            data = response.json()
            return data["PropertyTable"]["Properties"][0]
        except (KeyError, IndexError, requests.exceptions.JSONDecodeError):
            st.warning("Unable to parse PubChem compound summary.")
            return None
    else:
        st.error(f"Error fetching compound summary from PubChem. Status code: {response.status_code}")
        return None

# Helper: Fetch PubChem Similar Compounds
def fetch_similar_compounds(cid, threshold=90):
    """
    Retrieve similar compounds for a given PubChem CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/similarity/cids/JSON"
    params = {"Threshold": threshold}
    response = requests.get(url, params=params)
    if response.status_code == 200:
        try:
            data = response.json()
            return data["IdentifierList"]["CID"]
        except KeyError:
            st.warning("No similar compounds found.")
            return []
    else:
        st.error(f"Error fetching similar compounds. Status code: {response.status_code}")
        return []

# Helper: Save results to Excel
def save_to_excel(data, filename="results.xlsx"):
    """
    Save data to an Excel file.
    """
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
    st.title("Comprehensive Compound Analysis")
    st.write("Analyze compounds using PubChem for bioactivity and target data.")

    # Input: SMILES string
    smiles = st.text_input("Enter the SMILES string of your compound:")

    if st.button("Analyze Compound"):
        if smiles:
            # Validate SMILES
            if not validate_smiles(smiles):
                st.error("Invalid SMILES string. Please provide a valid structure.")
                return

            # Fetch CID from PubChem
            with st.spinner("Fetching PubChem CID..."):
                cid = get_cid_from_smiles(smiles)

            if cid:
                st.success(f"CID retrieved: {cid}")

                # Fetch PubChem compound summary
                with st.spinner("Fetching compound summary from PubChem..."):
                    compound_summary = fetch_pubchem_summary(cid)
                    if compound_summary:
                        st.subheader("PubChem Compound Summary")
                        st.json(compound_summary)

                # Fetch similar compounds from PubChem
                with st.spinner("Fetching similar compounds from PubChem..."):
                    similar_compounds = fetch_similar_compounds(cid)
                    if similar_compounds:
                        st.subheader("Similar Compounds (PubChem)")
                        st.write(similar_compounds)
                    else:
                        st.warning("No similar compounds found.")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()

