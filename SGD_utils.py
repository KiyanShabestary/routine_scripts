import pandas as pd
import json
import requests
from datetime import datetime


ORGANISM_ID = "559292"

def get_gene_description(SGD_gene_ID):
	'''
	Get single gene description using API. 

	Input:
		- gene_SGD_ID <str> : SGD gene ID to retrieve information from ("SGD:SXXXXXXXXX").

	Output:
		- data <dict> : Gene description for specified gene.
	'''

	# Define the base URL and endpoint
	base_url = "https://www.alliancegenome.org"
	endpoint = "/api/gene/{id}"

	# Construct the full URL by replacing {id} in the endpoint with the actual gene ID
	url = f"{base_url}{endpoint}".replace("{id}", SGD_gene_ID)

	# Make the GET request to the API
	response = requests.get(url)

	# Check if the request was successful
	if response.status_code == 200:
		# Parse and print the JSON response
		original_data = response.json()
		
		# Decompose and retain essential data only
		simplified_data = {
			"id": original_data.get("id"),
			"symbol": original_data.get("symbol"),
			"name": original_data.get("name"),
			"synonym1": original_data.get("synonyms")[0] if len(original_data.get("synonyms", [])) > 0 else None,
			"synonym2": original_data.get("synonyms")[1] if len(original_data.get("synonyms", [])) > 1 else None,
			"synonym3": original_data.get("synonyms")[2] if len(original_data.get("synonyms", [])) > 2 else None,
			"chromosome": original_data.get("genomeLocations")[0].get("chromosome") if original_data.get("genomeLocations") else None,
			"start": original_data.get("genomeLocations")[0].get("start") if original_data.get("genomeLocations") else None,
			"end": original_data.get("genomeLocations")[0].get("end") if original_data.get("genomeLocations") else None,
			"strand": original_data.get("genomeLocations")[0].get("strand") if original_data.get("genomeLocations") else None,
			"uniprotID": None,
			"geneDescription": original_data.get("geneSynopsis"),
			"geneDescriptionAuto": original_data.get("automatedGeneSynopsis")
		}

		# Extract Uniprot ID from crossReferenceMap's 'other' list
		for ref in original_data.get("crossReferenceMap", {}).get("other", []):
			if "UniProtKB" in ref.get("name", ""):
				simplified_data["uniprotID"] = ref.get("name").split(":")[-1]
				break

		return simplified_data

	elif response.status_code == 404:
		print(f"No data found for gene ID: {SGD_gene_ID}")
		return 0
	else:
		# Print error message if the request was not successful
		print(f"Failed to retrieve data: {response.status_code}")
		print(response.text)
		return 0


def get_gene_table(taxon_ID=ORGANISM_ID):
	'''
	Get gene IDs for specific taxon_ID. 

	Input:
		- taxon_ID <num> : taxon ID (NNNNNN).

	Output:
		- data <pandas> : Gene table of all genes for specified taxon ID.
	'''

	# Define the base URL and endpoint
	base_url = "https://www.alliancegenome.org"
	endpoint = "/api/geneMap"

	# Construct the full URL by replacing {id} in the endpoint with the actual gene ID
	url = f"{base_url}{endpoint}"

	# Define the query parameters
	params = {
		"taxonID": taxon_ID,  
		"page": 1,  # Default page number
		"rows": 10000  # Set to arbitrary high number to retrieve all entries
	}

	# Make the GET request to the API
	#response = requests.get(url, headers=headers)
	response = requests.get(url, params=params)

	# Check if the request was successful
	if response.status_code == 200:
		# Parse and print the JSON response
		data = response.json()

		entries_dict = {}

		for gene in data["results"]:
			entries_dict[gene["id"]] = get_gene_description(gene["id"])

		gene_table = pd.DataFrame.from_dict(entries_dict, orient='index')

		return gene_table
	else:
		# Print error message if the request was not successful
		print(f"Failed to retrieve data: {response.status_code}")
		print(response.text)
		return 0


def export_gene_table(taxon_ID=559292):
	'''
	Export gene table in a csv to work offline and more efficiently retrieve gene description.
	'''

	retrieve_date = datetime.now().date()

	path = f"SGD_S288C_gene_table_{retrieve_date}.csv"
	table = get_gene_table(taxon_ID=taxon_ID)

	table.to_csv(path, index=False)

	return 0

def import_gene_table(path):
	'''
	Import gene table from local file since calling the API and assembling at each iteration can be time-consuming.
	'''
	return pd.read_csv(path)
