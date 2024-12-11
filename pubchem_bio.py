import requests
import pandas as pd
import streamlit as st
from rdkit import Chem
from io import BytesIO
import urllib.parse

# Helper: Validate SMILES
def validate_smiles(smiles):
    """
    Basic validation for SMILES string.
    """
    return bool(smiles.strip())


# Helper: Fetch PubChem CID
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

# Helper: Fetch PubChem BioAssay Data
def fetch_pubchem_bioassay(cid):
    """
    Retrieve bioassay data for a given PubChem CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if "AssaySummary" in data:
            assays = data["AssaySummary"]["Assay"]
            results = [
                {
                    "Assay ID": assay.get("AID"),
                    "Description": assay.get("Description", "N/A"),
                    "Outcome": assay.get("Outcome", "N/A"),
                    "Target Name": assay.get("TargetName", "N/A"),
                }
                for assay in assays
            ]
            return results
    st.warning("No bioassay data found in PubChem.")
    return []

# Helper: Fetch ChEMBL Targets
def fetch_chembl_targets(smiles):
    """
    Fetch targets for a compound from ChEMBL using SMILES.
    """
    smiles_encoded = urllib.parse.quote(smiles)  # Encode SMILES for URL
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?smiles={smiles_encoded}"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            data = response.json()
            results = [
                {
                    "Molecule ChEMBL ID": molecule["molecule_chembl_id"],
                    "Molecule Name": molecule.get("pref_name", "N/A"),
                }
                for molecule in data.get("molecules", [])
            ]
            return results
        except requests.exceptions.JSONDecodeError:
            st.warning("Received a non-JSON response from ChEMBL.")
            st.write("Raw Response:", response.text)
            return []
    else:
        st.error(f"Error fetching targets from ChEMBL. Status code: {response.status_code}")
        st.write("Response content:", response.text)
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
    st.title("Compound Bioactivity and Target Finder")
    st.write("Analyze compounds using PubChem and ChEMBL for bioactivity data and targets.")

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

                # Fetch Canonical SMILES
                with st.spinner("Fetching canonical SMILES from PubChem..."):
                    canonical_smiles = fetch_canonical_smiles(cid)
                    if canonical_smiles:
                        st.write(f"Canonical SMILES from PubChem: {canonical_smiles}")
                        smiles = canonical_smiles  # Use canonical SMILES for downstream queries

                # Fetch BioAssay data from PubChem
                with st.spinner("Fetching bioassay data from PubChem..."):
                    pubchem_bioassay = fetch_pubchem_bioassay(cid)
                    if pubchem_bioassay:
                        st.subheader("BioAssay Data (PubChem)")
                        st.dataframe(pubchem_bioassay)

                        # Save PubChem bioassay data to Excel
                        excel_file = save_to_excel(pubchem_bioassay, "pubchem_bioassay.xlsx")
                        if excel_file:
                            st.download_button(
                                label="Download PubChem BioAssay Data as Excel",
                                data=excel_file,
                                file_name="pubchem_bioassay.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )
                    else:
                        st.warning("No bioassay data found in PubChem for this compound.")

                # Fetch protein targets from ChEMBL
                with st.spinner("Fetching targets from ChEMBL..."):
                    chembl_targets = fetch_chembl_targets(smiles)
                    if chembl_targets:
                        st.subheader("Protein Targets (ChEMBL)")
                        st.dataframe(chembl_targets)

                        # Save ChEMBL target data to Excel
                        chembl_excel = save_to_excel(chembl_targets, "chembl_targets.xlsx")
                        if chembl_excel:
                            st.download_button(
                                label="Download ChEMBL Targets as Excel",
                                data=chembl_excel,
                                file_name="chembl_targets.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )
                    else:
                        st.warning("No targets found in ChEMBL for this compound.")
                        st.write("You can manually verify the compound on ChEMBL:", "https://www.ebi.ac.uk/chembl/")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")


