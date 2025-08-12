import os
import sqlite3
import json
import requests
import urllib.parse
from gui_app.normalize_compound import get_compound_data

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATABASE_PATH = os.path.join(BASE_DIR, 'database.db')
SCHEMA_PATH = os.path.join(BASE_DIR, 'schema.sql')

# List of common antibiotics for seeding
ANTIBIOTICS_TO_ADD = [
    "Penicillin G", "Amoxicillin", "Tetracycline", "Ciprofloxacin",
    "Azithromycin", "Doxycycline", "Metronidazole", "Clindamycin",
    "Vancomycin", "Streptomycin", "Erythromycin", "Cephalexin",
    "Sulfamethoxazole", "Trimethoprim", "Levofloxacin", "Gentamicin",
    "Ampicillin", "Chloramphenicol", "Rifampicin", "Isoniazid",
    "Cefazolin", "Cefoxitin", "Cefuroxime", "Cefotaxime", "Ceftriaxone",
    "Cefepime", "Ceftaroline", "Imipenem", "Meropenem", "Ertapenem",
    "Doripenem", "Aztreonam", "Polymyxin B", "Colistin", "Bacitracin",
    "Neomycin", "Kanamycin", "Tobramycin", "Amikacin", "Netilmicin",
    "Spectinomycin", "Linezolid", "Tedizolid", "Daptomycin", "Fidaxomicin",
    "Rifaximin", "Telavancin", "Dalbavancin", "Oritavancin", "Fosfomycin",
    "Nitrofurantoin", "Methenamine", "Phenazopyridine", "Mupirocin", "Retapamulin",
    "Fusidic acid", "Gramicidin", "Tyrocidine", "Valinomycin", "Nisin",
    "Polymyxin E", "Capreomycin", "Cycloserine", "Ethambutol", "Pyrazinamide",
    "Dapsone", "Clofazimine", "Bedaquiline", "Delamanid", "Pretomanid",
    "Thiacetazone", "Ethionamide", "Prothionamide", "Para-aminosalicylic acid",
    "Viomycin", "Kanamycin A", "Amikacin", "Arbekacin", "Dibekacin",
    "Sisomicin", "Isepamicin", "Plazomicin", "Apramycin", "Hygromycin B",
    "Puromycin", "Streptothricin", "Bleomycin", "Actinomycin D", "Mithramycin",
    "Doxorubicin", "Daunorubicin", "Epirubicin", "Idarubicin", "Mitoxantrone",
    "Bleomycin A2", "Bleomycin B2", "Plicamycin", "Valrubicin", "Zorubicin"
]

def get_db_connection():
    """Establishes a connection to the SQLite database."""
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        return conn
    except sqlite3.Error as e:
        print(f"Database connection error: {e}")
        return None

def init_db():
    """Initializes the database from the schema.sql file."""
    print(f"Initializing database at {DATABASE_PATH}...")
    if os.path.exists(DATABASE_PATH):
        os.remove(DATABASE_PATH)
        print("Removed existing database file.")

    conn = get_db_connection()
    if not conn:
        print("Could not create database connection. Aborting.")
        return

    try:
        with open(SCHEMA_PATH, 'r') as f:
            schema_sql = f.read()
        conn.executescript(schema_sql)
        conn.commit()
        print("Database schema created successfully.")
    except sqlite3.Error as e:
        print(f"Error creating schema: {e}")
    except FileNotFoundError:
        print(f"Schema file not found at {SCHEMA_PATH}")
    finally:
        conn.close()

def fetch_compounds_by_name_list(compound_names: list, limit: int = 100) -> list:
    """Fetches compounds by a list of names and filters for ChEMBL/PubChem sources."""
    compounds = []
    print(f"Searching for compounds from provided list...")
    for name in compound_names:
        if len(compounds) >= limit:
            break
        print(f"Processing: {name}")
        compound_data = get_compound_data(name)

        if 'error' not in compound_data and compound_data.get('inchi_key'):
            compounds.append(compound_data)
            print(f"  -> Added {compound_data.get('common_name', compound_data.get('iupac_name', 'Unknown'))} (InChIKey found).")
        else:
            error_msg = compound_data.get('error', 'No InChIKey found')
            print(f"  -> Skipped {name} (Error: {error_msg}).")

    print(f"Finished fetching compounds. Total found: {len(compounds)}")
    return compounds

def seed_database():
    """Connects to the database and inserts the compound list."""
    print("\nStarting database seeding process...")
    conn = get_db_connection()
    if not conn:
        print("Database connection failed. Aborting.")
        return

    inserted_count = 0
    skipped_count = 0

    try:
        cursor = conn.cursor()

        # Fetch compounds dynamically
        compounds_to_insert = fetch_compounds_by_name_list(ANTIBIOTICS_TO_ADD, limit=100)

        for data in compounds_to_insert:
            print(f"Processing: {data.get('common_name', data.get('iupac_name', 'Unknown'))}")

            # Check if the compound already exists using inchi_key
            cursor.execute("SELECT compound_id FROM compounds WHERE inchi_key = ?", (data['inchi_key'],))
            if cursor.fetchone():
                print("  -> Compound already exists. Skipping.")
                skipped_count += 1
                continue

            # Prepare sources for storage
            sources_json = json.dumps(data.get('sources', []))

            # Insert the new compound
            sql = """
            INSERT INTO compounds (
                iupac_name, common_name, smiles_raw, smiles_normalized, inchi, inchi_key, 
                molecular_formula, molecular_weight, fingerprint, structure_2d_svg, metadata
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            val = (
                data.get('iupac_name'),
                data.get('common_name'),
                data['smiles_raw'], 
                data['smiles_normalized'], 
                data['inchi'], 
                data['inchi_key'],
                data['molecular_formula'], 
                data['molecular_weight'], 
                data['fingerprint'], 
                data['structure_2d_svg'],
                sources_json # Store sources as JSON in metadata
            )
            cursor.execute(sql, val)
            inserted_count += 1
            print("  -> Successfully inserted.")

        conn.commit()

    except sqlite3.Error as err:
        print(f"A database error occurred: {err}")
        conn.rollback()
    finally:
        conn.close()
        print("\nSeeding process finished.")
        print(f"Inserted {inserted_count} new compounds.")
        print(f"Skipped {skipped_count} existing compounds.")

if __name__ == "__main__":
    init_db() # Create the database and schema first
    seed_database() # Then populate it with data
