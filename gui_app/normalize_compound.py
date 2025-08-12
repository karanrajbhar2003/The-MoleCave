import argparse
import requests
import urllib.parse
import re
from typing import Union
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import rdMolDescriptors
from molvs import standardize_smiles
import urllib3
import json
import logging
from rdkit.Chem.Draw import MolDraw2DSVG # Import MolDraw2DSVG

# Configure logging for this module
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Simple input detectors
INCHIKEY_REGEX = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

def looks_like_inchikey(text: str) -> bool:
    return bool(INCHIKEY_REGEX.match(text.strip()))

def looks_like_pubchem_cid(text: str) -> bool:
    return text.isdigit()

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

def get_inchikey_from_pubchem(cid: str) -> Union[str, None]:
    try:
        response = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON",
            verify=False
        )
        if response.status_code == 200:
            data = response.json()
            properties = data['PropertyTable']['Properties'][0]
            return properties.get('InChIKey')
        else:
            logger.debug(f"PubChem InChIKey API failed for CID {cid}: {response.status_code} - {response.text}")
    except (requests.exceptions.RequestException, KeyError) as e:
        logger.error(f"Exception in get_inchikey_from_pubchem for CID {cid}: {e}")
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

def generate_2d_structure_svg(smiles: str) -> Union[str, None]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            drawer = MolDraw2DSVG(300, 300)
            drawer.drawOptions().addStereoAnnotation = True
            drawer.drawOptions().addAtomIndices = False
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            return drawer.GetDrawingText()
    except Exception as e:
        logger.error(f"Error generating 2D structure SVG for SMILES {smiles}: {e}")
    return None

# --- UniChem API Functions for DrugBank mapping ---

_unichem_drugbank_src_id_cache: Union[str, None] = None

def _get_unichem_drugbank_src_id() -> Union[str, None]:
    global _unichem_drugbank_src_id_cache
    if _unichem_drugbank_src_id_cache is not None:
        logger.debug("UniChem DrugBank source ID from cache.")
        return _unichem_drugbank_src_id_cache
    try:
        resp = requests.get("https://www.ebi.ac.uk/unichem/rest/sources", timeout=15, verify=False)
        if resp.status_code != 200:
            logger.debug(f"UniChem sources API failed: {resp.status_code} - {resp.text}")
            return None
        sources = resp.json()
        # Find DrugBank source id by name match
        for src in sources:
            name = str(src.get('name', '')).lower()
            if 'drugbank' in name:
                _unichem_drugbank_src_id_cache = str(src.get('src_id'))
                logger.debug(f"UniChem DrugBank source ID found: {_unichem_drugbank_src_id_cache}")
                break
        if not _unichem_drugbank_src_id_cache:
            logger.debug("DrugBank source not found in UniChem sources list")
        return _unichem_drugbank_src_id_cache
    except requests.exceptions.RequestException as e:
        logger.error(f"Exception fetching UniChem sources: {e}")
        return None


def get_drugbank_url_from_inchikey(inchikey: str) -> Union[str, None]:
    logger.debug(f"Attempting to get DrugBank URL for InChIKey: {inchikey}")
    try:
        src_id = _get_unichem_drugbank_src_id()
        logger.debug(f"UniChem DrugBank source ID used: {src_id}")
        if not src_id:
            logger.debug("No UniChem DrugBank source ID available.")
            return None
        # Query UniChem for all sources mapped to this InChIKey
        resp = requests.get(f"https://www.ebi.ac.uk/unichem/rest/inchikey/{inchikey}", timeout=15, verify=False)
        logger.debug(f"UniChem inchikey API response status for {inchikey}: {resp.status_code}")
        if resp.status_code != 200:
            logger.debug(f"UniChem inchikey API failed for {inchikey}: {resp.status_code} - {resp.text}")
            return None
        mappings = resp.json() or []
        logger.debug(f"UniChem mappings for {inchikey}: {json.dumps(mappings, indent=2)}")
        # mappings expected: list of {'src_id': '2', 'src_compound_id': 'DB01050', ...}
        db_id = None
        for m in mappings:
            if str(m.get('src_id')) == str(src_id):
                db_id = m.get('src_compound_id')
                if db_id:
                    logger.debug(f"DrugBank ID found in mappings: {db_id}")
                    break
        if not db_id:
            logger.debug(f"No DrugBank ID found in UniChem mappings for {inchikey}.")
            return None
        # Build modern DrugBank URL
        drugbank_url = f"https://go.drugbank.com/drugs/{db_id}"
        logger.debug(f"Constructed DrugBank URL: {drugbank_url}")
        return drugbank_url
    except (requests.exceptions.RequestException, ValueError) as e:
        logger.error(f"Exception in get_drugbank_url_from_inchikey for {inchikey}: {e}")
        return None

def get_drugbank_url_from_pubchem_cid(cid: str) -> Union[str, None]:
    """Fallback: Use PubChem xrefs to get DrugBank accession from a CID."""
    # Try dedicated DrugBank xrefs endpoint if supported
    try:
        resp = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/DrugBank/JSON",
            verify=False,
            timeout=20,
        )
        if resp.status_code == 200:
            data = resp.json()
            # Expected structure: InformationList -> Information[0] -> DrugBank -> ["DBxxxx"]
            info_list = data.get('InformationList', {}).get('Information', [])
            if info_list:
                db_ids = info_list[0].get('DrugBank') or []
                if db_ids:
                    return f"https://go.drugbank.com/drugs/{db_ids[0]}"
    except (requests.exceptions.RequestException, ValueError, KeyError) as e:
        logger.debug(f"PubChem DrugBank xrefs endpoint failed for CID {cid}: {e}")

    # Fallback to RegistryID filtered by DrugBank source parameter
    try:
        resp = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON?source=DrugBank",
            verify=False,
            timeout=20,
        )
        if resp.status_code == 200:
            data = resp.json()
            info_list = data.get('InformationList', {}).get('Information', [])
            if info_list:
                ids = info_list[0].get('RegistryID') or []
                # Typically contains the DB accession
                for reg in ids:
                    if isinstance(reg, str) and reg.upper().startswith('DB'):
                        return f"https://go.drugbank.com/drugs/{reg}"
    except (requests.exceptions.RequestException, ValueError, KeyError) as e:
        logger.debug(f"PubChem RegistryID xrefs fallback failed for CID {cid}: {e}")

    return None

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
    inchi = None
    inchi_key = None
    structure_2d_svg = None # Re-initialize
    molecular_formula = None # Re-initialize
    molecular_weight = None # Re-initialize
    smiles_normalized = None # Re-initialize

    # Helper to add a source only once
    def add_source(db_name: str, url: str) -> None:
        nonlocal sources
        if not any(s.get('db_name') == db_name for s in sources):
            sources.append({'db_name': db_name, 'url': url})

    # 1. Try to parse identifier as SMILES
    try:
        mol_from_smiles = Chem.MolFromSmiles(identifier)
        if mol_from_smiles:
            smiles = identifier
            add_source('User-provided', '#')
            logger.debug(f"Identifier recognized as SMILES: {smiles}")
    except Exception as e:
        logger.debug(f"Failed to parse identifier as SMILES: {e}")

    # 2. Attempt PubChem search (by name or InChIKey if SMILES was found)
    current_pubchem_cid = None
    # If identifier looks like a PubChem CID, use it directly
    if looks_like_pubchem_cid(identifier):
        current_pubchem_cid = identifier

    # If identifier looks like an InChIKey, set it and resolve CID directly
    if not current_pubchem_cid and looks_like_inchikey(identifier):
        inchi_key = identifier
        try:
            current_pubchem_cid = get_cid_from_inchikey(inchi_key)
        except Exception as e:
            logger.debug(f"Failed to resolve CID from provided InChIKey {inchi_key}: {e}")

    if smiles: # If SMILES was provided, try to get CID from its InChIKey
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                inchi = Chem.MolToInchi(mol)
                inchi_key = Chem.InchiToInchiKey(inchi)
                logger.debug(f"Derived InChIKey from SMILES: {inchi_key}")
                current_pubchem_cid = get_cid_from_inchikey(inchi_key)
        except Exception as e:
            logger.error(f"Error getting InChIKey/CID from SMILES: {e}")
    
    if not current_pubchem_cid: # If no CID yet, try by identifier name
        current_pubchem_cid = get_cid_from_pubchem(identifier)

    if current_pubchem_cid:
        logger.debug(f"PubChem CID {current_pubchem_cid} found for {identifier}.")
        pubchem_cid = current_pubchem_cid # Set the main pubchem_cid
        # Always add PubChem source link when CID is available
        add_source('PubChem', f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_cid}")

        # Ensure we fetch InChIKey from PubChem whenever CID is available
        if not inchi_key:
            fetched_inchi_key = get_inchikey_from_pubchem(pubchem_cid)
            if fetched_inchi_key:
                inchi_key = fetched_inchi_key

        # Try to get SMILES and names from PubChem if CID is available
        pubchem_smiles = get_smiles_from_pubchem(pubchem_cid)
        if pubchem_smiles:
            if not smiles: # Only set smiles if not already found
                smiles = pubchem_smiles
            
            # Get names from PubChem
            names = get_names_from_pubchem(pubchem_cid)
            if names['iupac_name']:
                iupac_name = names['iupac_name']
            if names['common_name']:
                common_name = names['common_name']
        else:
            logger.debug(f"No SMILES found from PubChem for CID {pubchem_cid}.")

    # 3. Attempt ChEMBL search
    chembl_smiles = get_smiles_from_chembl(identifier)
    if chembl_smiles:
        if not smiles: # Only set smiles if not already found
            smiles = chembl_smiles
        add_source('ChEMBL', f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={urllib.parse.quote(identifier)}")

    # 4. After all SMILES attempts, if SMILES is available, generate SVG and derive InChI/InChIKey
    if smiles:
        structure_2d_svg = generate_2d_structure_svg(smiles)
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                inchi = Chem.MolToInchi(mol)
                inchi_key = Chem.InchiToInchiKey(inchi)
                logger.debug(f"Derived InChIKey: {inchi_key}")
        except Exception as e:
            logger.error(f"Error deriving InChI/InChIKey from SMILES {smiles}: {e}")

    # 5. Populate DrugBank source if possible via UniChem
    logger.debug(f"Attempting DrugBank URL resolution with InChIKey: {inchi_key}")
    if inchi_key:
        drugbank_url = get_drugbank_url_from_inchikey(inchi_key)
        if drugbank_url:
            add_source('DrugBank', drugbank_url)
            logger.debug(f"DrugBank URL successfully added: {drugbank_url}")
        else:
            logger.debug("DrugBank URL not found via UniChem for this InChIKey.")
    else:
        logger.debug("InChIKey not available for DrugBank URL resolution.")

    # 5b. If UniChem did not yield DrugBank, try PubChem xrefs fallback
    if not any(s['db_name'] == 'DrugBank' for s in sources) and pubchem_cid:
        db_url_from_pubchem = get_drugbank_url_from_pubchem_cid(pubchem_cid)
        if db_url_from_pubchem:
            add_source('DrugBank', db_url_from_pubchem)
            logger.debug(f"DrugBank URL successfully added via PubChem xrefs: {db_url_from_pubchem}")

    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
                molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
                logger.debug(f"Calculated Molecular Formula: {molecular_formula}, Molecular Weight: {molecular_weight}")
        except Exception as e:
            logger.error(f"Error calculating molecular formula or weight for SMILES {smiles}: {e}")

    smiles_normalized = None
    fingerprint = None # Initialize fingerprint
    if smiles:
        try:
            smiles_normalized = standardize_smiles(smiles)
            mol = Chem.MolFromSmiles(smiles_normalized) # Use normalized SMILES for fingerprint
            if mol:
                # Calculate Morgan fingerprint (ECFP4, 2048 bits)
                fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048).ToBitString()
        except Exception as e:
            logger.error(f"Error standardizing SMILES or calculating fingerprint for {smiles}: {e}")

    return {
        "smiles_raw": smiles,
        "smiles_normalized": smiles_normalized,
        "iupac_name": iupac_name,
        "common_name": common_name,
        "sources": sources,
        "pubchem_cid": pubchem_cid,
        "structure_2d_svg": structure_2d_svg,
        "inchi": inchi,
        "inchi_key": inchi_key,
        "molecular_formula": molecular_formula,
        "molecular_weight": molecular_weight,
        "fingerprint": fingerprint # Add fingerprint to the returned dictionary
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch and process chemical data.")
    parser.add_argument("identifier", type=str, help="The chemical identifier (name, SMILES, etc.).")
    args = parser.parse_args()

    result = get_compound_data(args.identifier)
    print(json.dumps(result, indent=4))
