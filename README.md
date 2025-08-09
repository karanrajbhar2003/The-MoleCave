
# Chemical Compound Repository

This project is a web-based application for storing, searching, and viewing chemical compound data. It uses Flask for the backend and a self-contained SQLite database, making it highly portable and easy to set up.

## Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution) must be installed. This is the recommended way to install the `rdkit` dependency, especially on Windows.

## Setup & Running the Application

Follow these steps to get the application running on a new machine.

### 1. Create the Conda Environment

Open your terminal or Anaconda Prompt and navigate to the project's root directory. Run the following command to create a new Conda environment using the provided file. This will install all required libraries.

```bash
conda env create -f environment.yml
```

### 2. Activate the Environment

Before running any scripts, you must activate the new environment.

```bash
conda activate chemical-repo
```

### 3. Create and Seed the Database

Run the following command from the project's root directory. This single command will:

1.  Create a new `database.db` file.
2.  Set up the necessary tables using `schema.sql`.
3.  Populate the database with a list of common chemical compounds.

```bash
python seed.py
```

### 4. Run the Web Application

Finally, start the Flask web server with this command:

```bash
python -m gui_app.app
```

## Usage

-   Once the server is running, open your web browser and go to `http://127.0.0.1:5000`.
-   You can now:
    -   **View** the list of pre-loaded compounds.
    -   **Search** for compounds by name, SMILES, or formula.
    -   **Upload** new compounds by providing their SMILES string.

