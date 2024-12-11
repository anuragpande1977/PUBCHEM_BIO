import requests
import pandas as pd
import streamlit as st
from io import BytesIO
from rdkit import Chem

# Helper: Validate SMILES
def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

# Helper: Fetch CID from PubChem
def get_cid_from_structure(smiles):
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
        st.error(f"Error fetching CID. Status code: {response.status_code}")
        return None

# Helper: Fetch bioassay data
def get_bioassay_data(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'AssaySummary' in data:
            return data['AssaySummary']
    return None

# Helper: Fetch predictions from PASS Online
def fetch_pass_predictions(smiles):
    url = "http://www.pharmaexpert.ru/passonline/api"
    payload = {"smiles": smiles}
    try:
        response = requests.post(url, json=payload, timeout=10)  # Add timeout
        response.raise_for_status()  # Raise an error for HTTP codes 4xx/5xx
        
        # Attempt to parse JSON
        try:
            return response.json()
        except requests.exceptions.JSONDecodeError:
            st.error("Error: The API response is not in JSON format.")
            st.write("Response content:", response.text)  # Debugging
            return None
    except requests.exceptions.RequestException as e:
        st.error(f"API request failed: {e}")
        return None

# Helper: Save results to Excel
def save_to_excel(data, filename="results.xlsx"):
    output = BytesIO()
    try:
        df = pd.DataFrame(data)
        df.to_excel(output, index=False, engine='openpyxl')  # Save to Excel buffer
        output.seek(0)  # Reset pointer to the beginning of the stream
        return output
    except Exception as e:
        st.error(f"Error saving to Excel: {e}")
        return None

# Streamlit main app
def main():
    st.title("Chemical Analysis: PubChem BioAssay and PASS Predictions")
    st.write("Analyze compounds using PubChem and PASS Online for bioassay activities and predicted biological activity.")

    # Initialize session state
    if "cid" not in st.session_state:
        st.session_state.cid = None
    if "bioassay_data" not in st.session_state:
        st.session_state.bioassay_data = None
    if "pass_results" not in st.session_state:
        st.session_state.pass_results = None

    # Input: SMILES string
    smiles = st.text_input("Enter the SMILES string of your compound:")

    if st.button("Analyze Compound"):
        if smiles:
            # Validate SMILES
            if not validate_smiles(smiles):
                st.error("Invalid SMILES string. Please enter a valid structure.")
                return
            
            # Fetch CID
            with st.spinner("Fetching PubChem CID..."):
                st.session_state.cid = get_cid_from_structure(smiles)

            if st.session_state.cid:
                st.success(f"CID retrieved: {st.session_state.cid}")

                # Fetch BioAssay Data
                with st.spinner("Fetching BioAssay data..."):
                    st.session_state.bioassay_data = get_bioassay_data(st.session_state.cid)

                # Fetch PASS Predictions
                with st.spinner("Fetching PASS Online predictions..."):
                    st.session_state.pass_results = fetch_pass_predictions(smiles)

                st.success("Analysis complete!")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

    # Display results
    if st.session_state.cid:
        # BioAssay Data
        st.subheader("PubChem BioAssay Data")
        if st.session_state.bioassay_data:
            st.json(st.session_state.bioassay_data)
        else:
            st.warning("No BioAssay data found for this compound.")

        # PASS Predictions
        st.subheader("PASS Online Predictions")
        if st.session_state.pass_results:
            pass_data = pd.DataFrame(st.session_state.pass_results)
            st.dataframe(pass_data)

            # Save PASS Predictions to Excel
            excel_file = save_to_excel(st.session_state.pass_results, "pass_predictions.xlsx")
            if excel_file:
                st.download_button(
                    label="Download PASS Predictions as Excel",
                    data=excel_file,
                    file_name="pass_predictions.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )
        else:
            st.warning("No PASS predictions available for this compound.")

if __name__ == "__main__":
    main()

