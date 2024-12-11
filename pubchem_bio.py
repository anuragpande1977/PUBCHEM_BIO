import requests
import pandas as pd
import streamlit as st

def get_cid_from_structure(smiles):
    """
    Convert a SMILES string to a PubChem CID with debug statements.
    """
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"
    params = {'smiles': smiles}
    response = requests.get(url, params=params)
    
    st.write(f"SMILES to CID API Response: {response.status_code}, {response.json()}")  # Debug
    
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

def get_similar_cids(cid, threshold):
    """
    Retrieve CIDs for compounds similar to the given CID with debug statements.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/similarity/cids/JSON"
    params = {'Threshold': threshold}
    response = requests.get(url, params=params)
    
    st.write(f"Similarity API Response for CID {cid}: {response.status_code}, {response.json()}")  # Debug
    
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID']
        else:
            st.warning("No similar compounds found for the given CID.")
            return []
    else:
        st.error(f"Error fetching similar compounds. Status code: {response.status_code}")
        return []

def get_bioassay_data(cid):
    """
    Retrieve bioassay data for a given PubChem CID with debug statements.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)
    
    st.write(f"BioAssay API Response for CID {cid}: {response.status_code}, {response.json()}")  # Debug
    
    if response.status_code == 200:
        data = response.json()
        if 'AssaySummary' in data:
            return data['AssaySummary']
        else:
            st.warning(f"No 'AssaySummary' key found in the response for CID {cid}.")
            return None
    else:
        st.error(f"Error fetching bioassay data. Status code: {response.status_code}")
        return None

def save_to_excel(results, filename="bioassay_results.xlsx"):
    """
    Save bioassay results to an Excel file with debug statements.
    """
    try:
        df = pd.DataFrame(results)
        df.to_excel(filename, index=False)
        st.write(f"Results saved to {filename}")  # Debug
        return filename
    except Exception as e:
        st.error(f"Error saving to Excel: {e}")
        return None

def main():
    """
    Streamlit app for PubChem BioAssay Finder with debugging enabled.
    """
    st.title("PubChem BioAssay Finder with Debugging")
    st.write("Find potential biological targets for your compound using PubChem's BioAssay database.")

    smiles = st.text_input("Enter the SMILES string of your compound:")
    threshold = st.slider("Set Similarity Threshold (0-100):", min_value=0, max_value=100, value=80)
    
    if st.button("Search"):
        if smiles:
            with st.spinner("Fetching CID from PubChem..."):
                cid = get_cid_from_structure(smiles)
            
            if cid:
                st.success(f"CID retrieved: {cid}")
                
                # Step 1: Fetch BioAssay data for exact CID
                with st.spinner("Fetching BioAssay data for exact CID..."):
                    bioassay_data = get_bioassay_data(cid)
                
                results = []
                
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
                        results.append({"CID": cid, "Assay": name, "Target": target, "Outcome": outcome})
                else:
                    st.warning("No bioassay data found for the exact CID.")
                
                # Step 2: Fetch similar CIDs
                with st.spinner(f"Fetching similar compounds with threshold {threshold}..."):
                    similar_cids = get_similar_cids(cid, threshold)
                
                if similar_cids:
                    st.write("### BioAssay Data for Similar Compounds:")
                    for similar_cid in similar_cids[:10]:  # Limit to top 10 similar compounds
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
                                results.append({"CID": similar_cid, "Assay": name, "Target": target, "Outcome": outcome})
                        else:
                            st.write(f"No bioassay data found for CID {similar_cid}.")
                else:
                    st.warning("No similar compounds found.")
                
                # Step 3: Save results to Excel
                if results:
                    st.write("### Save Results")
                    filename = save_to_excel(results)
                    if filename:
                        with open(filename, "rb") as file:
                            st.download_button(label="Download Results as Excel", data=file, file_name="bioassay_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()

