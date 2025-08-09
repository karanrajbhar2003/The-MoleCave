-- Chemical Repository Database Schema for SQLite

-- Table for storing core chemical information
CREATE TABLE IF NOT EXISTS compounds (
    compound_id INTEGER PRIMARY KEY AUTOINCREMENT,
    iupac_name TEXT,
    common_name TEXT,
    molecular_formula TEXT,
    molecular_weight REAL,
    smiles_raw TEXT,
    smiles_normalized TEXT,
    inchi TEXT,
    inchi_key TEXT UNIQUE,
    source_url TEXT,
    structure_2d_svg TEXT,
    structure_3d_pdb TEXT,
    fingerprint TEXT,
    source_db TEXT,
    metadata TEXT, -- Storing JSON as TEXT in SQLite
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Indexes for efficient searching
CREATE INDEX IF NOT EXISTS idx_inchi_key ON compounds(inchi_key);
CREATE INDEX IF NOT EXISTS idx_smiles_normalized ON compounds(smiles_normalized);
