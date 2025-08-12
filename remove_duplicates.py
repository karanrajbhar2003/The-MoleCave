import sqlite3

def remove_duplicates():
    conn = sqlite3.connect('database.db')
    cursor = conn.cursor()

    # Find all InChI Keys that have more than one entry
    cursor.execute("""
        SELECT inchi_key
        FROM compounds
        GROUP BY inchi_key
        HAVING COUNT(inchi_key) > 1
    """)
    duplicate_keys = cursor.fetchall()

    if not duplicate_keys:
        print("No duplicate compounds found.")
        conn.close()
        return

    print(f"Found {len(duplicate_keys)} duplicate InChI Keys. Cleaning database...")

    for key_tuple in duplicate_keys:
        inchi_key = key_tuple[0]
        print(f"  Processing duplicate InChI Key: {inchi_key}")

        # Get all compound_ids for the duplicate InChI Key, ordered by ID
        cursor.execute("""
            SELECT compound_id
            FROM compounds
            WHERE inchi_key = ?
            ORDER BY compound_id
        """, (inchi_key,))
        compound_ids = [row[0] for row in cursor.fetchall()]

        # Keep the first one, delete the rest
        ids_to_delete = compound_ids[1:]
        
        if ids_to_delete:
            print(f"    Keeping compound_id {compound_ids[0]} and deleting IDs: {ids_to_delete}")
            cursor.execute(f"""
                DELETE FROM compounds
                WHERE compound_id IN ({','.join('?' for _ in ids_to_delete)})
            """, ids_to_delete)
            conn.commit()
        else:
            print(f"    Found duplicate key but no extra records to delete for {inchi_key}. This is unusual.")

    print("\nDatabase cleaning complete.")
    conn.close()

if __name__ == '__main__':
    remove_duplicates()