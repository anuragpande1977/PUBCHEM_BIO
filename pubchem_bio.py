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
    
    st.write(f"Similarity API Response: {response.status_code}, {response.json()}")  # Debugging output
    
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
    Streamlit app for PubChem BioAssay Finder with session state management.
    """
    st.title("PubChem BioAssay Finder")
    st.write("Find potential biological targets for your compound using PubChem's BioAssay database.")
    
    # Initialize session state
    if "exact_results" not in st.session_state:
        st.session_state.exact_results = None
    if "exact_found" not in st.session_state:
        st.session_state.exact_found = False
    if "cid" not in st.session_state:
        st.session_state.cid = None
    if "similar_results" not in st.session_state:
        st.session_state.similar_results = None
    
    # Input SMILES
    smiles = st.text_input("Enter the SMILES string of your compound:")
    
    # Exact search button
    if st.button("Search for Exact BioAssay Data"):
        if smiles:
            with st.spinner("Fetching CID from PubChem..."):
                st.session_state.cid = get_cid_from_structure(smiles)
            
            if st.session_state.cid:
                st.success(f"CID retrieved: {st.session_state.cid}")
                with st.spinner("Fetching BioAssay data for exact CID..."):
                    bioassay_data = get_bioassay_data(st.session_state.cid)
                    
                    if bioassay_data:
                        st.session_state.exact_found = True
                        if bioassay_data['type'] == 'AssaySummary':
                            st.session_state.exact_results = bioassay_data['data'].get('Assay', [])
                        elif bioassay_data['type'] == 'Table':
                            st.session_state.exact_results = bioassay_data['data']
                    else:
                        st.session_state.exact_found = False
                        st.warning("No bioassay data found for the exact CID.")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")
    
    # Show exact results if found
    if st.session_state.exact_found:
        activity_count = len(st.session_state.exact_results)
        st.write(f"### Exact BioAssay Activities Found: {activity_count}")
        filename = save_to_excel(st.session_state.exact_results, "exact_bioassay_results.xlsx")
        if filename:
            with open(filename, "rb") as file:
                st.download_button(label="Download Exact Results as Excel", data=file, file_name=filename, mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    else:
        # Similarity threshold search
        threshold = st.slider("Set Similarity Threshold (0-100):", min_value=0, max_value=100, value=80, key="threshold")
        if st.button("Search for Similar Compounds"):
            if st.session_state.cid:
                with st.spinner(f"Fetching similar compounds with threshold {threshold}..."):
                    similar_cids = get_similar_cids(st.session_state.cid, threshold)
                    if similar_cids:
                        st.session_state.similar_results = []
                        for similar_cid in similar_cids[:10]:  # Limit to top 10 similar compounds
                            bioassay_data = get_bioassay_data(similar_cid)
                            if bioassay_data:
                                if bioassay_data['type'] == 'AssaySummary':
                                    st.session_state.similar_results.extend(bioassay_data['data'].get('Assay', []))
                                elif bioassay_data['type'] == 'Table':
                                    st.session_state.similar_results.extend(bioassay_data['data'])
                        st.write(f"### Total BioAssay Activities Found for Similar Compounds: {len(st.session_state.similar_results)}")
                        filename = save_to_excel(st.session_state.similar_results, "similar_bioassay_results.xlsx")
                        if filename:
                            with open(filename, "rb") as file:
                                st.download_button(label="Download Similar Results as Excel", data=file, file_name=filename, mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                    else:
                        st.warning("No similar compounds found.")
            else:
                st.warning("No CID available. Perform an exact search first.")

if __name__ == "__main__":
    main()
