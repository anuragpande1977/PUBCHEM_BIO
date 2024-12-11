import requests
import pandas as pd
import streamlit as st
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

# Helper: Fetch NCBI Target Genes/Proteins
def fetch_ncbi_gene_protein(cid):
    """
    Retrieve associated genes and proteins using NCBI Entrez utilities.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/targets/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            data = response.json()
            results = [
                {"Gene/Protein": target["TargetName"], "Target Type": target.get("TargetType", "N/A")}
                for target in data.get("Targets", [])
            ]
            return results
        except KeyError:
            st.warning("No gene or protein targets found.")
            return []
    else:
        st.error(f"Error fetching target data from NCBI. Status code: {response.status_code}")
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
    st.title("Compound Analysis with Activity Predictions and Targets")
    st.write("Analyze compounds using PubChem and NCBI for activity and target data.")

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

                # Fetch PubChem BioAssay Data
                with st.spinner("Fetching bioassay data from PubChem..."):
                    bioassay_data = fetch_pubchem_bioassay(cid)
                    if bioassay_data:
                        st.subheader("BioAssay Data (PubChem)")
                        st.dataframe(bioassay_data)

                        # Save to Excel
                        excel_file = save_to_excel(bioassay_data, "bioassay_data.xlsx")
                        if excel_file:
                            st.download_button(
                                label="Download BioAssay Data as Excel",
                                data=excel_file,
                                file_name="bioassay_data.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )

                # Fetch NCBI Target Data
                with st.spinner("Fetching target genes and proteins from NCBI..."):
                    target_data = fetch_ncbi_gene_protein(cid)
                    if target_data:
                        st.subheader("Target Genes/Proteins (NCBI)")
                        st.dataframe(target_data)

                        # Save to Excel
                        excel_file = save_to_excel(target_data, "target_data.xlsx")
                        if excel_file:
                            st.download_button(
                                label="Download Target Data as Excel",
                                data=excel_file,
                                file_name="target_data.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()


