from flask import Flask, render_template, request, redirect, url_for, flash, jsonify
import sqlite3
import os
import json

# Import the new function from the local module
from .normalize_compound import get_compound_data

app = Flask(__name__)
# A secret key is needed for flashing messages
app.secret_key = 'a_random_secret_key'

# --- Database Configuration ---
# Get the absolute path for the database file
DATABASE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'database.db')

def get_db_connection():
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        # This allows accessing columns by name
        conn.row_factory = sqlite3.Row
        return conn
    except sqlite3.Error as e:
        flash(f"Database connection error: {e}", "danger")
        return None

# --- Routes ---

@app.route('/')
def index():
    
    conn = get_db_connection()
    if not conn:
        return render_template('index.html', count=0, error="Database connection failed.")

    count = 0
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM compounds")
        # fetchone() for sqlite returns a tuple, access by index
        result = cursor.fetchone()
        if result:
            count = result[0]
    except sqlite3.Error as err:
        flash(f"Error fetching data: {err}")
        count = 0
    finally:
        if conn:
            conn.close()

    return render_template('index.html', count=count)

@app.route('/upload', methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':
        identifier = request.form.get('identifier', '').strip()

        if not identifier:
            flash("Compound identifier is required.", "warning")
            return render_template('upload.html')

        # Process the identifier (this now handles names, SMILES, etc.)
        try:
            data = get_compound_data(identifier)
        except Exception as e:
            flash(f"Error processing compound data: {e}", "danger")
            print(f"DEBUG: Exception in get_compound_data: {e}")
            return render_template('upload.html')

        if 'error' in data and not data.get('smiles_raw'):
            if data.get('sources'):
                flash(f"Could not fully process compound: {data['error']}. However, external links were found and saved.", "info")
            else:
                flash(f"Error processing identifier: {data['error']}", "danger")
                print(f"DEBUG: Error processing identifier: {data['error']}")
                return render_template('upload.html')

        # If no smiles_raw is present, we cannot proceed with full compound data, but we might have sources.
        if not data.get('smiles_raw') and not data.get('sources'):
            flash(f"Could not process '{identifier}' and no external sources were found.", "danger")
            print(f"DEBUG: No SMILES or sources found for '{identifier}'")
            return render_template('upload.html')

        conn = get_db_connection()
        if not conn:
            flash("Database connection failed.", "danger")
            return render_template('upload.html')

        try:
            cursor = conn.cursor()
            existing_compound = None

            # 1. Check for duplicates using InChIKey (most reliable)
            if data.get('inchi_key'):
                cursor.execute("SELECT compound_id FROM compounds WHERE inchi_key = ?", (data['inchi_key'],))
                existing_compound = cursor.fetchone()
                if existing_compound:
                    flash(f"This compound (InChIKey match) already exists in the database (ID: {existing_compound['compound_id']}).", "info")
                    return redirect(url_for('compound_details', compound_id=existing_compound['compound_id']))

            # 2. If no InChIKey match, check for duplicates using smiles_normalized
            if not existing_compound and data.get('smiles_normalized'):
                cursor.execute("SELECT compound_id FROM compounds WHERE smiles_normalized = ?", (data['smiles_normalized'],))
                existing_compound = cursor.fetchone()
                if existing_compound:
                    flash(f"This compound (normalized SMILES match) already exists in the database (ID: {existing_compound['compound_id']}).", "info")
                    return redirect(url_for('compound_details', compound_id=existing_compound['compound_id']))

            # 3. If no normalized SMILES match, check for duplicates using smiles_raw
            if not existing_compound and data.get('smiles_raw'):
                cursor.execute("SELECT compound_id FROM compounds WHERE smiles_raw = ?", (data['smiles_raw'],))
                existing_compound = cursor.fetchone()
                if existing_compound:
                    flash(f"This compound (raw SMILES match) already exists in the database (ID: {existing_compound['compound_id']}).", "info")
                    return redirect(url_for('compound_details', compound_id=existing_compound['compound_id']))

            # 4. If no other match, check for duplicates using common_name
            if not existing_compound and data.get('common_name'):
                cursor.execute("SELECT compound_id FROM compounds WHERE common_name = ? COLLATE NOCASE", (data['common_name'],))
                existing_compound = cursor.fetchone()
                if existing_compound:
                    flash(f"This compound (common name match) already exists in the database (ID: {existing_compound['compound_id']}).", "info")
                    return redirect(url_for('compound_details', compound_id=existing_compound['compound_id']))

            # Prepare the SQL INSERT statement
            sql = """
            INSERT INTO compounds (
                iupac_name, common_name, smiles_raw, smiles_normalized, inchi, inchi_key,
                molecular_formula, molecular_weight, fingerprint, structure_2d_svg, metadata
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            val = (
                data.get('iupac_name'),
                data.get('common_name'),
                data.get('smiles_raw'), # Use .get() as it might be missing if only sources found
                data.get('smiles_normalized'), # Use .get() as it might be missing
                data.get('inchi'), # Use .get() as it might be missing
                data.get('inchi_key'), # Use .get() as it might be missing
                data.get('molecular_formula'), # Use .get() as it might be missing
                data.get('molecular_weight'), # Use .get() as it might be missing
                data.get('fingerprint'), # Use .get() as it might be missing
                data.get('structure_2d_svg'), # Use .get() as it might be missing
                json.dumps(data.get('sources', []))
            )
            cursor.execute(sql, val)
            conn.commit()
            flash("Compound added successfully!", "success")
        except sqlite3.Error as err:
            flash(f"Database error: {err}", "danger")
        finally:
            if conn:
                conn.close()

        # If a compound was added, redirect to its details page if possible
        if data.get('inchi_key'):
            conn = get_db_connection()
            if conn:
                cursor = conn.cursor()
                cursor.execute("SELECT compound_id FROM compounds WHERE inchi_key = ?", (data['inchi_key'],))
                new_compound = cursor.fetchone()
                if new_compound:
                    return redirect(url_for('compound_details', compound_id=new_compound['compound_id']))
            flash("Compound added, but could not retrieve details page.", "info")
        return redirect(url_for('index'))

    return render_template('upload.html')

@app.route('/search')
def search():
    conn = get_db_connection()
    if not conn:
        flash("Database connection failed.")
        return render_template('search.html', compounds=[])

    compounds = []
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT compound_id, iupac_name, common_name, smiles_normalized, molecular_formula FROM compounds ORDER BY common_name COLLATE NOCASE ASC, iupac_name COLLATE NOCASE ASC")
        compounds = cursor.fetchall()
    except sqlite3.Error as err:
        flash(f"Error fetching compounds: {err}")
    finally:
        if conn:
            conn.close()

    return render_template('search.html', compounds=compounds)

@app.route('/search_api')
def search_api():
    query = request.args.get('query', '').strip()
    conn = get_db_connection()
    if not conn:
        return jsonify(error="Database connection failed."), 500

    results = []
    try:
        cursor = conn.cursor()
        if query:
            like_query = f"%{query.lower()}%"
            sql = """
            SELECT compound_id, iupac_name, common_name, smiles_normalized, molecular_formula
            FROM compounds
            WHERE LOWER(iupac_name) LIKE ?
              OR LOWER(common_name) LIKE ?
              OR smiles_normalized LIKE ?
              OR smiles_raw LIKE ?
              OR inchi_key LIKE ?
              OR molecular_formula LIKE ?
            ORDER BY common_name COLLATE NOCASE ASC, iupac_name COLLATE NOCASE ASC
            """
            cursor.execute(sql, (like_query, like_query, like_query, like_query, like_query, like_query))
        else:
            cursor.execute("SELECT compound_id, iupac_name, common_name, smiles_normalized, molecular_formula FROM compounds ORDER BY common_name COLLATE NOCASE ASC, iupac_name COLLATE NOCASE ASC")
        
        results = [dict(row) for row in cursor.fetchall()]

    except sqlite3.Error as err:
        return jsonify(error=f"Search error: {err}"), 500
    finally:
        if conn:
            conn.close()

    return jsonify(results)

@app.route('/compound/<int:compound_id>')
def compound_details(compound_id):
    print(f"DEBUG: Entering compound_details route for compound_id: {compound_id}")
    conn = get_db_connection()
    if not conn:
        print("DEBUG: Database connection failed in compound_details.")
        return redirect(url_for('index'))

    compound = None
    sources = []
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM compounds WHERE compound_id = ?", (compound_id,))
        compound = cursor.fetchone()
        if not compound:
            flash("Compound not found.")
            print(f"DEBUG: Compound with ID {compound_id} not found in database.")
            return redirect(url_for('index'))
        
        # Convert sqlite3.Row to a dictionary for easier JSON serialization
        compound_dict = dict(compound)

        if compound_dict['metadata']:
            try:
                sources = json.loads(compound_dict['metadata'])
                print(f"DEBUG: Successfully parsed metadata for compound {compound_id}.")
            except json.JSONDecodeError as e:
                flash("Error parsing metadata for compound.", "warning")
                print(f"DEBUG: Error parsing metadata for compound {compound_id}: {e}")

    except sqlite3.Error as err:
        flash(f"Error fetching compound details: {err}")
        print(f"DEBUG: SQLite error fetching compound details for {compound_id}: {err}")
        return redirect(url_for('index'))
    finally:
        if conn:
            conn.close()

    print(f"DEBUG: Rendering compound.html for compound {compound_id}. Compound data: {compound_dict}")
    return render_template('compound.html', compound=compound_dict, sources=sources)

@app.route('/about')
def about():
    return render_template('about.html')

if __name__ == '__main__':
    # Use 0.0.0.0 to make it accessible from the network
    app.run(host='0.0.0.0', port=5000, debug=True)