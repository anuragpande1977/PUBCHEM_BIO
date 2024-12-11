import requests
import pandas as pd
import streamlit as st
from io import BytesIO

# Helper: Fetch similar compounds from ZINC20
def fetch_similar_compounds_from_zinc(smiles, threshold=0.8):
    url = f"https://zinc20.docking.org/substances.txt?structure={smiles}&threshold={threshold}"
    response = requests.get(url)
    if response.status_code == 200:
        compounds = response.text.splitlines()
        return compounds[:10]  # Return top 10 similar compounds
    else:
        st.error(f"Error fetching similar compounds from ZINC20. Status code: {response.status_code}")
        return []

# Helper: Fetch compound details from ZINC20
def fetch_compound_details(zinc_id):
    url = f"https://zinc20.docking.org/substances/{zinc_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        st.error(f"Error fetching details for ZINC ID {zinc_id}. Status code: {response.status_code}")
        return None

# Helper: Fetch protein targets from ZINC20
def fetch_targets_from_zinc(smiles):
    url = f"https://zinc20.docking.org/proteins.txt?structure={smiles}"
    response = requests.get(url)
    if response.status_code == 200:
        targets = response.text.splitlines()  # List of target protein names
        return targets
    else:
        st.error(f"Error fetching targets from ZINC20. Status code: {response.status_code}")
        return []

# Helper: Save results to Excel
def save_to_excel(data, filename="results.xlsx"):
    output = BytesIO()
    try:
        df = pd.DataFrame(data)
        df.to_excel(output, index=False, engine='openpyxl')
        output.seek(0)
        return output
    except Exception as e:
        st.error(f"Error saving to Excel: {e}")
        return None

# Streamlit main app
def main():
    st.title("ZINC20 Compound Similarity and Target Finder")
    st.write("Retrieve similar compounds, targets, and activity data using ZINC20.")

    # Input: SMILES string
    smiles = st.text_input("Enter the SMILES string of your compound:")

    if st.button("Analyze Compound"):
        if smiles:
            # Fetch similar compounds
            with st.spinner("Fetching similar compounds from ZINC20..."):
                similar_compounds = fetch_similar_compounds_from_zinc(smiles)

            # Display similar compounds
            if similar_compounds:
                st.subheader("Similar Compounds")
                st.write(similar_compounds)

                # Fetch and display compound details
                compound_details = []
                for zinc_id in similar_compounds:
                    details = fetch_compound_details(zinc_id)
                    if details:
                        compound_details.append(details)
                if compound_details:
                    st.subheader("Compound Details")
                    st.dataframe(compound_details)

                    # Save compound details to Excel
                    excel_file = save_to_excel(compound_details, "compound_details.xlsx")
                    if excel_file:
                        st.download_button(
                            label="Download Compound Details as Excel",
                            data=excel_file,
                            file_name="compound_details.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        )

            # Fetch protein targets
            with st.spinner("Fetching protein targets from ZINC20..."):
                targets = fetch_targets_from_zinc(smiles)
                if targets:
                    st.subheader("Protein Targets")
                    st.write(targets)
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()


