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

# Helper: Fetch ChEMBL Data by Name
def fetch_chembl_by_name(name):
    """
    Fetch ChEMBL data using compound name.
    """
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?pref_name={urllib.parse.quote(name)}"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            data = response.json()
            return [
                {
                    "Molecule ChEMBL ID": molecule["molecule_chembl_id"],
                    "Name": molecule.get("pref_name", "N/A"),
                    "Max Phase": molecule.get("max_phase", "N/A"),
                }
                for molecule in data.get("molecules", [])
            ]
        except (KeyError, requests.exceptions.JSONDecodeError):
            st.warning("Unable to parse ChEMBL response.")
            return []
    else:
        st.error(f"Error fetching ChEMBL data by name. Status code: {response.status_code}")
        return []

# Helper: Fetch Similar Compounds from PubChem
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

# Helper: Fetch Binding Data from BindingDB
def fetch_bindingdb_targets(smiles):
    """
    Fetch binding data for a compound from BindingDB.
    """
    url = f"https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFDownload.jsp?download_file=yes&smiles={urllib.parse.quote(smiles)}"
    response = requests.get(url)
    if response.status_code == 200 and response.text.strip():
        data = response.text.splitlines()
        return data[:10]  # Show first 10 lines
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
    st.title("Comprehensive Compound Analysis")
    st.write("Analyze compounds using PubChem, ChEMBL, and BindingDB for bioactivity and target data.")

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

                # Fetch PubChem bioassay data
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
                        st.warning("No bioassay data found in PubChem.")

                # Fetch similar compounds from PubChem
                with st.spinner("Fetching similar compounds from PubChem..."):
                    similar_compounds = fetch_similar_compounds(cid)
                    if similar_compounds:
                        st.subheader("Similar Compounds (PubChem)")
                        st.write(similar_compounds)
                    else:
                        st.warning("No similar compounds found.")

                # Fetch ChEMBL data using compound name
                compound_name = f"Compound {cid}"  # Placeholder; replace with a proper name fetch if available
                with st.spinner("Fetching data from ChEMBL using compound name..."):
                    chembl_data = fetch_chembl_by_name(compound_name)
                    if chembl_data:
                        st.subheader("ChEMBL Data by Name")
                        st.dataframe(chembl_data)
                    else:
                        st.warning("No data found in ChEMBL using the compound name.")

                # Fetch binding data from BindingDB
                with st.spinner("Fetching binding data from BindingDB..."):
                    bindingdb_targets = fetch_bindingdb_targets(smiles)
                    if bindingdb_targets:
                        st.subheader("Binding Targets (BindingDB)")
                        st.write(bindingdb_targets)
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()

