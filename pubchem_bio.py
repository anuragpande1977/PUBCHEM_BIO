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

# Helper: Fetch PASS Online Predictions
def fetch_pass_predictions(smiles):
    """
    Fetch predicted activities (PI and PA) using PASS Online.
    """
    url = "https://www.way2drug.com/PASSOnline/predict"
    params = {"smiles": smiles}
    response = requests.post(url, data=params)
    if response.status_code == 200:
        try:
            data = response.json()
            predictions = [
                {"Activity": activity, "PA": pred["PA"], "PI": pred["PI"]}
                for activity, pred in data["results"].items()
            ]
            return predictions
        except KeyError:
            st.warning("Unable to parse PASS Online predictions.")
            return []
    else:
        st.error(f"Error fetching predictions from PASS Online. Status code: {response.status_code}")
        return []

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

# Helper: Fetch PubChem Similar Compounds (Using SMILES)
def fetch_similar_compounds_by_smiles(smiles, threshold=80):
    """
    Retrieve similar compounds for a given SMILES string using PubChem.
    """
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/cids/JSON"
    params = {"smiles": smiles, "Threshold": threshold}
    response = requests.get(url, params=params)
    if response.status_code == 200:
        try:
            data = response.json()
            return data["IdentifierList"]["CID"]
        except KeyError:
            st.warning("No similar compounds found.")
            return []
    else:
        st.error(f"Error fetching similar compounds by SMILES. Status code: {response.status_code}")
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
    st.title("Comprehensive Compound Analysis with PI and PA Predictions")
    st.write("Analyze compounds using PubChem, PASS Online, and more.")

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

                # Fetch PASS Online predictions
                with st.spinner("Fetching activity predictions from PASS Online..."):
                    pass_predictions = fetch_pass_predictions(smiles)
                    if pass_predictions:
                        st.subheader("Predicted Activities (PASS Online)")
                        st.dataframe(pass_predictions)

                        # Save predictions to Excel
                        excel_file = save_to_excel(pass_predictions, "pass_predictions.xlsx")
                        if excel_file:
                            st.download_button(
                                label="Download PASS Predictions as Excel",
                                data=excel_file,
                                file_name="pass_predictions.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )
                    else:
                        st.warning("No predictions available from PASS Online.")

                # Fetch similar compounds from PubChem
                with st.spinner("Fetching similar compounds from PubChem using SMILES..."):
                    similar_compounds = fetch_similar_compounds_by_smiles(smiles)
                    if similar_compounds:
                        st.subheader("Similar Compounds (PubChem)")
                        st.write(similar_compounds)

                        # Save similar compounds to Excel
                        similar_data = [{"CID": cid} for cid in similar_compounds]
                        excel_file = save_to_excel(similar_data, "similar_compounds.xlsx")
                        if excel_file:
                            st.download_button(
                                label="Download Similar Compounds as Excel",
                                data=excel_file,
                                file_name="similar_compounds.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )
                    else:
                        st.warning("No similar compounds found.")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()


