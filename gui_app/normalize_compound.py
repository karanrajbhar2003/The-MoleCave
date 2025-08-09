import argparse
import requests
import urllib.parse
from typing import Union
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import rdMolDescriptors
from molvs import standardize_smiles
import urllib3
import json
import logging

# Configure logging for this module
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# --- PubChem API Functions ---

def get_cid_from_pubchem(identifier: str) -> Union[str, None]:
    try:
        encoded_identifier = urllib.parse.quote(identifier)
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_identifier}/cids/JSON",
            verify=False
        )
        if response.status_code == 200:
            data = response.json()
            return data['IdentifierList']['CID'][0]
        else:
            logger.debug(f"PubChem CID API failed for {identifier}: {response.status_code} - {response.text}")
    except (requests.exceptions.RequestException, KeyError) as e:
        logger.error(f"Exception in get_cid_from_pubchem for {identifier}: {e}")
    return None

def get_cid_from_inchikey(inchikey: str) -> Union[str, None]:
    try:
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON",
            verify=False
        )
        if response.status_code == 200:
            data = response.json()
            return data['IdentifierList']['CID'][0]
        else:
            logger.debug(f"PubChem CID from InChIKey API failed for {inchikey}: {response.status_code} - {response.text}")
    except (requests.exceptions.RequestException, KeyError) as e:
        logger.error(f"Exception in get_cid_from_inchikey for {inchikey}: {e}")
    return None

def get_smiles_from_pubchem(cid: str) -> Union[str, None]:
    try:
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,ConnectivitySMILES/JSON",
            verify=False
        )
        if response.status_code == 200:
            data = response.json()
            properties = data['PropertyTable']['Properties'][0]
            if 'CanonicalSMILES' in properties:
                return properties['CanonicalSMILES']
            elif 'ConnectivitySMILES' in properties:
                return properties['ConnectivitySMILES']
        else:
            logger.debug(f"PubChem SMILES API failed for CID {cid}: {response.status_code} - {response.text}")
    except (requests.exceptions.RequestException, KeyError) as e:
        logger.error(f"Exception in get_smiles_from_pubchem for CID {cid}: {e}")
    return None

def get_names_from_pubchem(cid: str) -> dict:
    iupac_name = None
    common_name = None
    try:
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON",
            verify=False
        )
        if response.status_code == 200:
            data = response.json()
            synonyms = data['InformationList']['Information'][0]['Synonym']
            if synonyms:
                common_name = synonyms[0]
        else:
            logger.debug(f"PubChem Synonyms API failed for CID {cid}: {response.status_code} - {response.text}")

        response_iupac = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON",
            verify=False
        )
        if response_iupac.status_code == 200:
            data_iupac = response_iupac.json()
            iupac_name = data_iupac['PropertyTable']['Properties'][0]['IUPACName']
        else:
            logger.debug(f"PubChem IUPAC API failed for CID {cid}: {response_iupac.status_code} - {response_iupac.text}")
    except (requests.exceptions.RequestException, KeyError) as e:
        logger.error(f"Exception in get_names_from_pubchem for CID {cid}: {e}")
    return {"iupac_name": iupac_name, "common_name": common_name}

# --- ChEMBL API Functions ---

def get_smiles_from_chembl(identifier: str) -> Union[str, None]:
    try:
        search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_synonyms__molecule_synonym__iexact={urllib.parse.quote(identifier)}"
        response = requests.get(search_url, verify=False)
        if response.status_code == 200:
            data = response.json()
            if data and 'molecules' in data and len(data['molecules']) > 0:
                molecule = data['molecules'][0]
                if (
                    'molecule_structures' in molecule
                    and molecule['molecule_structures'] is not None
                    and 'canonical_smiles' in molecule['molecule_structures']
                ):
                    return molecule['molecule_structures']['canonical_smiles']
            else:
                logger.debug(f"ChEMBL API returned no molecules for {identifier}")
        else:
            logger.debug(f"ChEMBL API failed for {identifier}: {response.status_code} - {response.text}")
    except (requests.exceptions.RequestException, KeyError) as e:
        logger.error(f"Exception in get_smiles_from_chembl for {identifier}: {e}")
    return None

# --- Core Data Generation Function ---

def get_compound_data(identifier: str) -> dict:
    logger.debug(f"Starting get_compound_data for identifier: {identifier}")
    if not identifier:
        logger.debug("Identifier is empty.")
        return {"error": "Identifier cannot be empty."}

    smiles = None
    iupac_name = None
    common_name = None
    sources = []
    pubchem_cid = None

    try:
        mol_from_smiles = Chem.MolFromSmiles(identifier)
        if mol_from_smiles:
            smiles = identifier
            sources.append({'db_name': 'User-provided', 'url': '#'})
            logger.debug(f"Identifier recognized as SMILES: {smiles}")
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    inchi = Chem.MolToInchi(mol)
                    inchi_key = Chem.InchiToInchiKey(inchi)
                    logger.debug(f"Derived InChIKey from SMILES: {inchi_key}")
                    pubchem_cid = get_cid_from_inchikey(inchi_key)
                    if pubchem_cid:
                        logger.debug(f"PubChem CID {pubchem_cid} found from InChIKey derived from SMILES.")
            except Exception as e:
                logger.error(f"Error getting InChIKey/CID from SMILES: {e}")
    except Exception as e:
        logger.error(f"Failed to parse identifier as SMILES: {e}")

    # Additional code would go here to handle name-based search and fallback
    return {
        "smiles": smiles,
        "iupac_name": iupac_name,
        "common_name": common_name,
        "sources": sources,
        "pubchem_cid": pubchem_cid
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch and process chemical data.")
    parser.add_argument("identifier", type=str, help="The chemical identifier (name, SMILES, etc.).")
    args = parser.parse_args()

    result = get_compound_data(args.identifier)
    logger.debug(json.dumps(result, indent=4))
