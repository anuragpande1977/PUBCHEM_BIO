import requests
import pandas as pd
import streamlit as st
from io import BytesIO

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

# Helper: Fetch bioassay data from PubChem
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

# Helper: Fetch protein targets from ChEMBL
def fetch_chembl_targets(smiles):
    """
    Fetch targets for a compound from ChEMBL using SMILES.
    """
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?smiles={smiles}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        targets = []
        for molecule in data.get("molecules", []):
            for activity in molecule.get("activities", []):
                if "target_chembl_id" in activity:
                    targets.append(activity["target_chembl_id"])
        return targets
    else:
        st.error(f"Error fetching targets from ChEMBL. Status code: {response.status_code}")
        return []

# Helper: Fetch binding data from BindingDB
def fetch_bindingdb_targets(smiles):
    """
    Fetch binding targets for a compound from BindingDB.
    """
    url = f"https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFDownload.jsp?download_file=yes&smiles={smiles}"
    response = requests.get(url)
    if response.status_code == 200:
        # BindingDB may return SDF or summary data; ensure proper processing
        data = response.text.splitlines()
        return data[:10]  # Show first 10 entries for brevity
    else:
        st.error(f"Error fetching binding data from BindingDB. Status code: {response.status_code}")
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

# Main Streamlit app
def main():
    st.title("Compound Bioactivity and Target Finder")
    st.write("Analyze compounds using PubChem, ChEMBL, and BindingDB for bioactivity data and targets.")

    # Input: SMILES string
    smiles = st.text_input("Enter the SMILES string of your compound:")

    if st.button("Analyze Compound"):
        if smiles:
            with st.spinner("Fetching PubChem CID..."):
                cid = get_cid_from_smiles(smiles)

            if cid:
                st.success(f"CID retrieved: {cid}")

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

                # Fetch protein targets from ChEMBL
                with st.spinner("Fetching targets from ChEMBL..."):
                    chembl_targets = fetch_chembl_targets(smiles)
                    if chembl_targets:
                        st.subheader("Protein Targets (ChEMBL)")
                        st.write(chembl_targets)

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

