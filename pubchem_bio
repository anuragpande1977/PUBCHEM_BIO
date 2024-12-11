import requests
import streamlit as st

# Function to get CID from SMILES
def get_cid_from_structure(smiles):
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
            return None
    else:
        return None

# Function to get BioAssay data
def get_bioassay_data(cid):
    """
    Retrieve bioassay data for a given PubChem CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        if 'AssaySummary' in data:
            return data['AssaySummary']
        else:
            return None
    else:
        return None

# Streamlit app
def main():
    st.title("PubChem BioAssay Finder")
    st.write("Find potential biological targets for your compound using PubChem's BioAssay database.")
    
    # User input
    smiles = st.text_input("Enter the SMILES string of your compound:")
    
    if st.button("Search"):
        if smiles:
            with st.spinner("Fetching CID from PubChem..."):
                cid = get_cid_from_structure(smiles)
            
            if cid:
                st.success(f"CID retrieved: {cid}")
                
                with st.spinner("Fetching BioAssay data..."):
                    bioassay_data = get_bioassay_data(cid)
                
                if bioassay_data:
                    st.write("### BioAssay Data Summary:")
                    for assay in bioassay_data.get('Assay', []):
                        name = assay.get('Name', 'Unknown')
                        target = assay.get('Target', {}).get('Description', 'Unknown target')
                        outcome = assay.get('Outcome', 'Unknown outcome')
                        st.write(f"- **Assay**: {name}")
                        st.write(f"  - **Target**: {target}")
                        st.write(f"  - **Outcome**: {outcome}")
                        st.write("---")
                else:
                    st.warning("No bioassay data found for this CID.")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()
