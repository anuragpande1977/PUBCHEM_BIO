import requests
import pandas as pd
import streamlit as st
from io import BytesIO

# Helper: Convert SMILES to PubChem CID
def get_cid_from_smiles(smiles):
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

# Helper: Retrieve BioAssay Data
def get_bioassay_data(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'AssaySummary' in data:
            return data['AssaySummary']['Assay']
    return None

# Helper: Fetch related PubMed articles
def fetch_pubmed_articles(query, max_results=10):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": max_results,
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        data = response.json()
        if 'esearchresult' in data and 'idlist' in data['esearchresult']:
            return data['esearchresult']['idlist']
    return []

# Helper: Retrieve PubMed article details
def fetch_pubmed_details(article_ids):
    if not article_ids:
        return []
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(article_ids),
        "retmode": "json",
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        data = response.json()
        if 'result' in data:
            results = data['result']
            return [
                {
                    "Title": results[article_id]["title"],
                    "Journal": results[article_id]["source"],
                    "Year": results[article_id]["pubdate"],
                }
                for article_id in article_ids if article_id in results
            ]
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
    st.title("Compound BioActivity Finder")
    st.write("Retrieve bioactivity data, targets, and related publications from NCBI.")

    # Input: SMILES string
    smiles = st.text_input("Enter the SMILES string of your compound:")

    if st.button("Analyze Compound"):
        if smiles:
            # Step 1: Fetch CID
            with st.spinner("Fetching PubChem CID..."):
                cid = get_cid_from_smiles(smiles)

            if cid:
                st.success(f"CID retrieved: {cid}")

                # Step 2: Fetch BioAssay Data
                with st.spinner("Fetching BioAssay data..."):
                    bioassay_data = get_bioassay_data(cid)

                # Step 3: Fetch Related Publications
                if bioassay_data:
                    st.subheader("BioAssay Data")
                    important_data = []
                    for assay in bioassay_data:
                        important_data.append({
                            "Assay ID": assay["AID"],
                            "Activity": assay.get("Outcome", "N/A"),
                            "Description": assay.get("Description", "N/A"),
                            "Target": assay.get("TargetName", "N/A")
                        })
                    st.dataframe(important_data)

                    # Save to Excel
                    excel_file = save_to_excel(important_data, "bioassay_data.xlsx")
                    if excel_file:
                        st.download_button(
                            label="Download BioAssay Data as Excel",
                            data=excel_file,
                            file_name="bioassay_data.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        )

                    # Fetch PubMed articles
                    st.spinner("Fetching related PubMed articles...")
                    targets = [assay["TargetName"] for assay in bioassay_data if "TargetName" in assay]
                    pubmed_articles = []
                    for target in targets:
                        articles = fetch_pubmed_articles(target)
                        pubmed_articles.extend(fetch_pubmed_details(articles))
                    
                    st.subheader("Related Publications")
                    st.dataframe(pubmed_articles)

                    # Save PubMed articles to Excel
                    pubmed_excel = save_to_excel(pubmed_articles, "pubmed_articles.xlsx")
                    if pubmed_excel:
                        st.download_button(
                            label="Download PubMed Articles as Excel",
                            data=pubmed_excel,
                            file_name="pubmed_articles.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        )
                else:
                    st.warning("No BioAssay data found for this compound.")
            else:
                st.error("Failed to retrieve CID. Please check your SMILES input.")
        else:
            st.error("Please enter a valid SMILES string.")

if __name__ == "__main__":
    main()


