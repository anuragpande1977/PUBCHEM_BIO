import requests
import streamlit as st

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

def get_similar_cids(cid):
    """
    Retrieve CIDs for compounds similar to the given CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/similarity/cids/JSON"
    params = {'Threshold': 90}  # Adjust similarity threshold (0-100)
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID']
        else:
            return []
    else:
        return []

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

def main():
    """
    Streamlit app for PubChem BioAssay Finder with structural similarity search.
    """
    st.title("PubChem BioAssay Finder")
    st.write("Find potential biological targets for your compound using PubChem's BioAssay database.")

    smiles = st.text_input("Enter the SMILES string of your compound:")

    if st.button("Search"):
        if smiles:
            with st.spinner("Fetching CID from PubChem..."):
                cid = get_cid_from_structure(smiles)
            
            if cid:
                st.success(f"CID retrieved: {cid}")
                
                # Step 1: Fetch BioAssay data for exact CID
                with st.spinner("Fetching BioAssay data for exact CID..."):
                    bioassay_data = get_bioassay_data(cid)
                
                if bioassay_data:
                    st.write("### BioAssay Data for Exact CID:")
                    for assay in bioassay_data.get('Assay', []):
                        name = assay.get('Name', 'Unknown')
                        target = assay.get('Target', {}).get('Description', 'Unknown target')
                        outcome = assay.get('Outcome', 'Unknown outcome')
                        st.write(f"- **Assay**: {name}")
                        st.write(f"  - **Target**: {target}")
                        st.write(f"  - **Outcome**: {outcome}")
                        st.write("---")
                else:
                    st.warning("No bioassay data found for the exact CID.")
                
                # Step 2: Fetch similar CIDs
                with st.spinner("Fetching similar compounds..."):
                    similar_cids = get_similar_cids(cid)
                
                if similar_cids:
                    st.write("### BioAssay Data for Similar Compounds:")
                    for similar_cid in similar_cids[:10]:  # Limit to top 10 similar compounds
                        st.write(f"Fetching bioassay data for similar CID: {similar_cid}")
                        bioassay_data = get_bioassay_data(similar_cid)
                        if bioassay_data:
                            for assay in bioassay_data.get('Assay', []):
                                name = assay.get('Name', 'Unknown')
                                target = assay.get('Target', {}).get('Description', 'Unknown target')
                                outcome = assay.get('Outcome', 'Unknown outcome')
                                st.write(f"- **Assay**: {name}")
                                st.write(f"  - **Target**: {target}")
                                st.write(f"  - **Outcome**: {outcome}")
                                st.write("---")
                        else:
                            st.write(f"No bioassay data found for CID {similar_cid}.")
                else:
                    st.warning("No similar compounds found.")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()

