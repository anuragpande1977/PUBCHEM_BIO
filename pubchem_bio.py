import requests
import pandas as pd
import streamlit as st
from io import BytesIO
from rdkit import Chem
import urllib.parse
response = requests.get(url)
if response.status_code == 200:
    try:
        data = response.json()
    except ValueError as e:
        st.error(f"JSON decode error: {e}")
        st.write("Raw Response Content:", response.text)  # Debugging info
        return []
else:
    st.error(f"Request failed. HTTP Status: {response.status_code}")
    st.write("Raw Response Content:", response.text)
    return []

def validate_and_canonicalize_smiles(smiles):
    """
    Validate and canonicalize the SMILES string using RDKit.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)  # Returns canonical SMILES
        else:
            st.error("Invalid SMILES string. Could not be parsed.")
            return None
    except Exception as e:
        st.error(f"Error during SMILES validation: {e}")
        return None

# Helper: Validate SMILES
def validate_smiles(smiles):
    """
    Basic validation for SMILES string. Ensure it's non-empty and well-formed.
    """
    return bool(smiles.strip())

# Helper: Canonicalize SMILES
def canonicalize_smiles(smiles):
    """
    Generate a canonical SMILES string using RDKit.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            st.error("Failed to canonicalize the SMILES string.")
            return None
    except Exception as e:
        st.error(f"Error canonicalizing SMILES: {e}")
        return None

# Helper: Fetch ChEMBL Molecule Info (Exact Match)
def fetch_chembl_exact(smiles):
    """
    Fetch molecule information from ChEMBL using the exact match endpoint.
    """
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?filter=molecule_structures__canonical_smiles__exact={urllib.parse.quote(smiles)}"
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
                st.warning("No molecules found in the ChEMBL exact match response.")
                return []
        except ValueError as e:
            st.error(f"JSON decode error: {e}")
            st.write("Raw Response Content:", response.text)
            return []
    else:
        st.error(f"Exact match failed. HTTP Status: {response.status_code}")
        st.write("Response Content:", response.text)
        return []

# Helper: Fetch ChEMBL Molecule Info (Substructure Search)
def fetch_chembl_substructure(smiles):
    """
    Fetch molecule information from ChEMBL using substructure search.
    """
    url = "https://www.ebi.ac.uk/chembl/api/data/substructure"
    headers = {"Content-Type": "application/json"}
    payload = {"smiles": smiles}  # Correct payload structure
    try:
        response = requests.post(url, headers=headers, json=payload)
        if response.status_code == 200:
            data = response.json()
            molecules = [
                {
                    "ChEMBL ID": molecule["molecule_chembl_id"],
                    "Name": molecule.get("pref_name", "N/A"),
                    "Max Phase": molecule.get("max_phase", "N/A"),
                }
                for molecule in data.get("molecules", [])
            ]
            if molecules:
                return molecules
            else:
                st.warning("No molecules found in the ChEMBL substructure response.")
                return []
        elif response.status_code == 400:
            st.error("Bad Request: Substructure search payload might be invalid.")
        else:
            st.error(f"Substructure search failed. HTTP Status: {response.status_code}")
            st.write("Response Content:", response.text)
    except Exception as e:
        st.error(f"Error during substructure search: {e}")
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

            # Canonicalize SMILES
            canonical_smiles = canonicalize_smiles(smiles)
            if not canonical_smiles:
                st.error("Failed to canonicalize SMILES. Analysis cannot proceed.")
                return

            # Fetch ChEMBL molecule information
            with st.spinner("Fetching molecule information from ChEMBL..."):
                chembl_data = fetch_chembl_exact(canonical_smiles)
                if chembl_data is None:  # Exact match failed, try substructure
                    chembl_data = fetch_chembl_substructure(canonical_smiles)
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
                bindingdb_targets = fetch_bindingdb_targets(canonical_smiles)
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



