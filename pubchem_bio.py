import requests
import pandas as pd
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
            st.warning("No CID found for the given SMILES string.")
            return None
    else:
        st.error(f"Error fetching CID. Status code: {response.status_code}")
        return None

def get_bioassay_data(cid):
    """
    Retrieve bioassay data for a given PubChem CID, handling multiple response formats.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        if 'AssaySummary' in data:
            return {'type': 'AssaySummary', 'data': data['AssaySummary']}
        if 'Table' in data:
            table = data['Table']
            columns = table.get('Columns', {}).get('Column', [])
            rows = table.get('Row', [])
            parsed_rows = [dict(zip(columns, row.get('Cell', []))) for row in rows]
            return {'type': 'Table', 'data': parsed_rows}
        return None
    else:
        return None

def get_similar_cids(cid, threshold):
    """
    Retrieve CIDs for compounds similar to the given CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/similarity/cids/JSON"
    params = {'Threshold': threshold}
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID']
        return []
    else:
        st.error(f"Error fetching similar compounds. Status code: {response.status_code}")
        return []

def save_to_excel(results, filename="bioassay_results.xlsx"):
    """
    Save bioassay results to an Excel file.
    """
    try:
        df = pd.DataFrame(results)
        df.to_excel(filename, index=False)
        return filename
    except Exception as e:
        st.error(f"Error saving to Excel: {e}")
        return None

def main():
    """
    Streamlit app for PubChem BioAssay Finder with fallback to similarity search.
    """
    st.title("PubChem BioAssay Finder")
    st.write("Find potential biological targets for your compound using PubChem's BioAssay database.")
    
    smiles = st.text_input("Enter the SMILES string of your compound:")
    
    if st.button("Search for Exact BioAssay Data"):
        if smiles:
            with st.spinner("Fetching CID from PubChem..."):
                cid = get_cid_from_structure(smiles)
            
            if cid:
                st.success(f"CID retrieved: {cid}")
                results = []
                activity_count = 0
                
                # Fetch BioAssay data for the exact compound
                with st.spinner("Fetching BioAssay data for exact CID..."):
                    bioassay_data = get_bioassay_data(cid)
                    if bioassay_data:
                        if bioassay_data['type'] == 'AssaySummary':
                            activity_count = len(bioassay_data['data'].get('Assay', []))
                            results = bioassay_data['data'].get('Assay', [])
                        elif bioassay_data['type'] == 'Table':
                            activity_count = len(bioassay_data['data'])
                            results = bioassay_data['data']
                
                if activity_count > 0:
                    st.write(f"### Exact BioAssay Activities Found: {activity_count}")
                    filename = save_to_excel(results)
                    if filename:
                        with open(filename, "rb") as file:
                            st.download_button(label="Download Exact Results as Excel", data=file, file_name="bioassay_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                else:
                    st.warning("No bioassay data found for the exact CID. Please adjust the similarity threshold to search for similar compounds.")
                    # Allow user to adjust threshold and search for similar compounds
                    threshold = st.slider("Set Similarity Threshold (0-100):", min_value=0, max_value=100, value=80)
                    if st.button("Search for Similar Compounds"):
                        with st.spinner(f"Fetching similar compounds with threshold {threshold}..."):
                            similar_cids = get_similar_cids(cid, threshold)
                            if similar_cids:
                                results = []
                                activity_count = 0
                                for similar_cid in similar_cids[:10]:  # Limit to top 10 similar compounds
                                    bioassay_data = get_bioassay_data(similar_cid)
                                    if bioassay_data:
                                        if bioassay_data['type'] == 'AssaySummary':
                                            activity_count += len(bioassay_data['data'].get('Assay', []))
                                            results.extend(bioassay_data['data'].get('Assay', []))
                                        elif bioassay_data['type'] == 'Table':
                                            activity_count += len(bioassay_data['data'])
                                            results.extend(bioassay_data['data'])
                                st.write(f"### Total BioAssay Activities Found for Similar Compounds: {activity_count}")
                                filename = save_to_excel(results)
                                if filename:
                                    with open(filename, "rb") as file:
                                        st.download_button(label="Download Similar Results as Excel", data=file, file_name="similar_bioassay_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                            else:
                                st.warning("No similar compounds found.")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()
