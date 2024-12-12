import requests
import pandas as pd
import streamlit as st
from io import BytesIO
import urllib.parse

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

# Helper: Fetch ChEMBL Molecule Info
def fetch_chembl_molecule_info(smiles):
    """
    Fetch molecule information from ChEMBL using the exact match endpoint.
    """
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search/{urllib.parse.quote(smiles)}"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            data = response.json()
            if "molecules" in data:
                molecules = [
                    {
                        "ChEMBL ID": molecule["molecule_chembl_id"],
                        "Name": molecule.get("pref_name", "N/A"),
                        "Max Phase": molecule.get("max_phase", "N/A"),
                    }
                    for molecule in data["molecules"]
                ]
                return molecules
            else:
                st.warning("No molecules found in the ChEMBL response.")
                return []
        except ValueError as e:
            st.error(f"JSON decode error: {e}")
            st.write("Raw Response Content:", response.text)  # Debugging info
            return []
    else:
        st.error(f"Failed to fetch data from ChEMBL. HTTP Status: {response.status_code}")
        st.write("Response Content:", response.text)  # Debugging info
        return []

# Helper: Fetch BindingDB Targets
def fetch_bindingdb_targets(smiles):
    """
    Fetch target binding data from BindingDB using SMILES.
    """
    url = f"https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFDownload.jsp?download_file=yes&smiles={urllib.parse.quote(smiles)}"
    response = requests.get(url)
    if response.status_code == 200 and response.text.strip():
        data = response.text.splitlines()
        return [{"Target Data": line} for line in data[:10]]  # Show first 10 entries
    else:
        st.warning("No binding data found in BindingDB for this compound.")
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
    st.title("Comprehensive Compound Analysis with ChEMBL and BindingDB Integration")
    st.write("Analyze compounds using PubChem, ChEMBL, and BindingDB.")

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

            # Fetch ChEMBL molecule information
            with st.spinner("Fetching molecule information from ChEMBL..."):
                chembl_data = fetch_chembl_molecule_info(smiles)
                if chembl_data:
                    st.subheader("ChEMBL Molecule Information")
                    st.dataframe(chembl_data)

                    # Save to Excel
                    excel_file = save_to_excel(chembl_data, "chembl_data.xlsx")
                    if excel_file:
                        st.download_button(
                            label="Download ChEMBL Data as Excel",
                            data=excel_file,
                            file_name="chembl_data.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        )
                else:
                    st.warning("No data found in ChEMBL.")

            # Fetch BindingDB target data
            with st.spinner("Fetching binding data from BindingDB..."):
                bindingdb_targets = fetch_bindingdb_targets(smiles)
                if bindingdb_targets:
                    st.subheader("Binding Targets (BindingDB)")
                    st.dataframe(bindingdb_targets)

                    # Save to Excel
                    excel_file = save_to_excel(bindingdb_targets, "bindingdb_targets.xlsx")
                    if excel_file:
                        st.download_button(
                            label="Download BindingDB Targets as Excel",
                            data=excel_file,
                            file_name="bindingdb_targets.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        )
                else:
                    st.warning("No binding data found in BindingDB.")

        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()



