import requests
import pandas as pd
import streamlit as st
import networkx as nx
import matplotlib.pyplot as plt

# Helper: Validate and process SMILES
def get_cid_from_structure(smiles):
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

# Helper: Retrieve bioassay data
def get_bioassay_data(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'AssaySummary' in data:
            return data['AssaySummary']
    return None

# Helper: Predict activity spectra (PASS Online)
def fetch_pass_predictions(smiles):
    url = "http://www.pharmaexpert.ru/passonline/api"
    payload = {"smiles": smiles}
    response = requests.post(url, json=payload)
    if response.status_code == 200:
        return response.json()
    else:
        st.error("Error fetching predictions from PASS Online.")
        return None

# Helper: SwissTargetPrediction for protein targets
def fetch_swiss_target_predictions(smiles):
    url = f"http://www.swisstargetprediction.ch/api/similarity/{smiles}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        st.error("Error fetching data from SwissTargetPrediction.")
        return None

# Helper: Visualize gene-protein interactions
def plot_gene_interactions(interactions):
    G = nx.Graph()
    for source, target in interactions:
        G.add_edge(source, target)

    plt.figure(figsize=(8, 6))
    nx.draw(G, with_labels=True, node_color="skyblue", font_weight="bold")
    st.pyplot(plt)

# Streamlit app setup
def main():
    st.title("PubChem BioAssay Finder and Target Prediction")
    st.write("Analyze compounds for bioassay activities, potential targets, and biological pathways.")

    # Initialize session state
    if "cid" not in st.session_state:
        st.session_state.cid = None
    if "bioassay_data" not in st.session_state:
        st.session_state.bioassay_data = None
    if "pass_results" not in st.session_state:
        st.session_state.pass_results = None
    if "targets" not in st.session_state:
        st.session_state.targets = None

    # Input SMILES string
    smiles = st.text_input("Enter the SMILES string of your compound:")
    
    if st.button("Analyze Compound"):
        if smiles:
            with st.spinner("Fetching PubChem CID..."):
                st.session_state.cid = get_cid_from_structure(smiles)

            if st.session_state.cid:
                st.success(f"CID retrieved: {st.session_state.cid}")
                with st.spinner("Fetching BioAssay data..."):
                    st.session_state.bioassay_data = get_bioassay_data(st.session_state.cid)

                with st.spinner("Predicting Activity Spectra (PASS Online)..."):
                    st.session_state.pass_results = fetch_pass_predictions(smiles)

                with st.spinner("Fetching SwissTargetPrediction data..."):
                    st.session_state.targets = fetch_swiss_target_predictions(smiles)

        else:
            st.error("Please enter a valid SMILES string.")
    
    # Display results
    if st.session_state.cid:
        st.subheader("PubChem BioAssay Data")
        if st.session_state.bioassay_data:
            st.json(st.session_state.bioassay_data)
        else:
            st.warning("No BioAssay data found for this compound.")

    if st.session_state.pass_results:
        st.subheader("Predicted Activity Spectra (PASS Online)")
        st.json(st.session_state.pass_results)

    if st.session_state.targets:
        st.subheader("Predicted Protein Targets (SwissTargetPrediction)")
        st.json(st.session_state.targets)

        # Example gene-protein interactions visualization
        st.subheader("Gene-Protein Interaction Map")
        example_interactions = [("Gene1", "ProteinA"), ("Gene2", "ProteinB"), ("Gene3", "ProteinC")]
        plot_gene_interactions(example_interactions)

if __name__ == "__main__":
    main()
