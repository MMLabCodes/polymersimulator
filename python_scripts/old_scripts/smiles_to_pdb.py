# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 14:30:07 2024

@author: danie
"""
import csv
import subprocess
import sys as sys
import itertools
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer
from rdkit.Chem.rdDistGeom import EmbedMolecule
from modules.sw_openmm import *
'''
  1. Here are functions defined to generate pdb files from smiles.
     The pdb files will be located in the above specified
         "pdb_file_destination"
     This also updates the residue_codes.csv to contain a specific residue code generated
         for each pdb.
    
    Note: There are a variety of functions defined here. However, all the functions defined are developed for use within the "SmilesToPDB" function. 
        It is recommend to only call the "SmilesToPDB" function. Although, there may be specific cases when the others are useful.
        
    Requirements:
        1. Some files are required for working with this function - namely .pdb files in the following path:
                    
            pdb_file_destination = current_working_directory + "/pdb_files/molecules"
            
        2. This python script located in the current_working_directory OR this python script executed in the current_working_directory by providing the correct path to it.

    IMPORTANT:
        Not currently in the mmlabcodes version of BCSW - but a different SMILES --> PDB should be available there.
        However, that version does not generate a unique residue code for each PDB file. The functions in this script are also imported into "csv-to_pdb.py" and utilised there.
'''
def SmilesToPDB_GenerateRescode(smiles, name, directory, residue_code_csv):
    """
    Converts a molecule specified by its SMILES representation to a PDB file.

    Parameters:
    - smiles (str): SMILES representation of the molecule.
    - name (str): Name of the molecule.
    - directory (str): Directory where the PDB file will be saved.
    - residue_code_csv (str): Path to the CSV file containing existing residue codes.

    Returns:
    None

    The function performs the following steps:
    1. Loads existing residue codes from the provided CSV file.
    2. Checks if the molecule name or SMILES is already in the database.
    3. If the entry exists, uses the existing residue code; otherwise, generates a unique 3-letter
       residue code excluding forbidden codes (which contains some examples "AAA", "BBB", "CCC" and amino acid residue codes.
    4. Updates the CSV file with the new entry if a new residue code is generated (If a new code isn't generated, the database already has info for that molecule).
    5. Converts the SMILES representation to a molecule object, adds hydrogens, and canonicalizes
       conformers.
    6. Replaces default "UNL" codes in the PDB content with the generated or existing residue code (if the molecule has already been assigned a residue code).
    7. Writes the PDB content to a file in the specified directory using the molecule's name.

    Note: The function utilizes various helper functions such as load_residue_codes,
    find_existing_entry, generate_unique_residue_code, update_residue_codes_csv, EmbedMolecule,
    and rdMolTransforms. These functions are defined and are available in this python file.
    """
    forbidden_codes = ["AAA", "BBB", "CCC", "UNL", "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "SEC", "TRP", "TYR", "VAL"]
    # Load existing residue codes from the CSV file
    residue_codes = load_residue_codes(residue_code_csv)
    # Check if the name or smiles is already in the database
    existing_entry = find_existing_entry(residue_codes, name, smiles)
    if existing_entry:
        # Use existing residue code
        residue_code = existing_entry[2]  # Assuming code is the third column
    else:
        # Generate a unique 3-letter residue code excluding forbidden codes
        residue_code = generate_unique_residue_code(residue_codes, forbidden_codes)
        # Update the CSV file with the new entry
        update_residue_codes_csv(name, smiles, residue_code, residue_code_csv)
    # Replace default "UNL" codes in the PDB content - NOTE, "UNL" is also a forbidden code as it is the default code.
    pdb_filepath = directory + "/" + name + ".pdb"
    smiles_to_pdb(smiles, pdb_filepath)    
    with open(pdb_filepath, "r") as pdb_file:
        lines = pdb_file.readlines()     
        for i, line in enumerate(lines):
            lines[i] = line.replace(" UNL", f" {residue_code}")   
    with open(pdb_filepath, "w") as pdb_file:              
        pdb_file.writelines(lines)
    return()

def load_residue_codes(residue_code_csv):
    """
    Loads existing residue codes from a CSV file.

    Parameters:
    - residue_code_csv (str): Path to the CSV file containing residue codes.

    Returns:
    list: A list of lists representing rows from the CSV file, where each row contains
          residue code information.

    The CSV file is expected to have rows containing residue code information.
    Each row may include data such as molecule name, SMILES representation, and
    the corresponding residue code. The function reads the CSV file and returns
    a list of lists, where each inner list represents a row from the CSV file.

    Example CSV structure:
    ```
    MoleculeName, SMILES, ResidueCode
    Example1, CCO, AAA
    Example2, CCN, BBB
    ...
    ```

    Note: Ensure that the CSV file has appropriate headers, and the function
    assumes that the first row contains column headers.
    """
    # Load existing residue codes from the CSV file
    residue_codes = []	
    with open(residue_code_csv, "r") as csv_file:
        reader = csv.reader(csv_file)
        residue_codes = [row for row in reader]
    return residue_codes

def generate_unique_residue_code(residue_codes, forbidden_codes=None):
    """
    Generates a unique 3-letter residue code not present in the database and not in the forbidden codes list.

    Parameters:
    - residue_codes (list): A list of lists representing existing residue code information.
                          Each inner list is expected to contain data such as molecule name,
                          SMILES representation, and the corresponding residue code.
    - forbidden_codes (list, optional): A list of 3-letter residue codes that are not allowed
                                      to be used. Defaults to None if no forbidden codes are specified..

    Returns:
    str: A unique 3-letter residue code.

    The function utilizes the existing residue codes and an optional list of forbidden codes
    to generate a new 3-letter residue code. It ensures that the generated code is not already
    present in the database and is not in the list of forbidden codes.

    If a unique code cannot be generated, a ValueError is raised.

    Note: The function assumes that the residue code is the third element in each inner list
    of the residue_codes parameter.
    """  
    # Generate a unique 3-letter residue code not already in the database and not in the forbidden codes list
    existing_codes = set(entry[2] for entry in residue_codes)  # Assuming code is the third column
    forbidden_codes = set(forbidden_codes) if forbidden_codes else set()

    # Generate all possible combinations of three letters
    all_combinations = [''.join(combination) for combination in itertools.product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=3)]

    # Find the first unused code not in the forbidden codes list
    new_code = next((code for code in all_combinations if code not in existing_codes and code not in forbidden_codes), None)

    if new_code is None:
        raise ValueError("Unable to generate a unique residue code.")

    return new_code

def update_residue_codes_csv(name, smiles, residue_code, residue_code_csv):
    """
    Updates a CSV file with a new entry if the entry does not already exist.

    Parameters:
    - name (str): Name of the molecule.
    - smiles (str): SMILES representation of the molecule.
    - residue_code (str): Residue code associated with the molecule.
    - residue_code_csv (str): Path to the CSV file containing existing residue codes.

    Returns:
    None

    The function checks if an entry with the provided name and smiles already exists
    in the CSV file. If not, it appends a new entry with the given information to the CSV file.

    Note: The function relies on helper functions such as load_residue_codes and find_existing_entry.
    Ensure these functions are defined and available in this python file.
    """
    # Update the CSV file with the new entry if it doesn't exist
    residue_codes = load_residue_codes(residue_code_csv)
    existing_entry = find_existing_entry(residue_codes, name, smiles)

    if not existing_entry:
        with open(residue_code_csv, "a", newline="") as csv_file:
            writer = csv.writer(csv_file)

            # Add a new entry to the CSV file
            writer.writerow([name, smiles, residue_code])

def find_existing_entry(residue_codes, name, smiles):
    """
    Finds an existing entry in a list of residue codes based on molecule name or SMILES representation.

    Parameters:
    - residue_codes (list): A list of lists representing existing residue code information.
                          Each inner list is expected to contain data such as molecule name,
                          SMILES representation, and the corresponding residue code.
    - name (str): Name of the molecule to search for.
    - smiles (str): SMILES representation of the molecule to search for.

    Returns:
    list or None: If an entry with the provided name or smiles is found, returns the corresponding
                 entry (a list). Otherwise, returns None.

    The function iterates through the list of residue codes and checks if any entry has a matching
    molecule name or SMILES representation. If a match is found, the corresponding entry is returned.
    If no match is found, None is returned.
    """
    for entry in residue_codes:
        if len(entry) >= 2 and (entry[0] == name or entry[1] == smiles):
            return entry
    return None

def pdb_to_mol(pdb_filename):
    """
    Converts a PDB file to a molecule object.

    Parameters:
    - pdb_filename (str): Path to the PDB file.

    Returns:
    Chem.Mol or None: A molecule object if conversion is successful, else None.

    The function reads the content of the PDB file, converts it to a molecule object,
    and assigns chiral tags based on the molecular structure. If the conversion is
    unsuccessful or the PDB file is empty, None is returned.
    """
    with open(pdb_filename, 'r') as pdb_file:
        pdb_content = pdb_file.read()

    mol = Chem.MolFromPDBBlock(pdb_content)
    if mol is not None:
        AllChem.AssignAtomChiralTagsFromStructure(mol)

    return mol

if __name__ == "__main__":
    '''
    Main script execution block:
    
    1. Check if <name> and <smiles> are provided as command-line arguments.
    
    NOTE: This will only run if this python script is executed - in other cases (see "csv_to_pdb.py") the functions are instead imported from this script
        and used there instead of this script being executed.
    '''
    if len(sys.argv) < 3:
        print("Error: Insufficient command-line arguments.")
        print("Usage: python3 script.py <name> <smiles>")
        print("Example: python3 smiles_to_pdb.py 'methane' 'C'")
        sys.exit(1)

    name = sys.argv[1]
    smiles = sys.argv[2]
    '''
    2. Get the current working directory as a variable called "current_working_directory".
       Note: This section utilizes bash commands executed with Python. To understand
             the behavior, run the command "pwd" in the Ubuntu or SCW shell.
    '''
    # Bash command to print the current working directory (pwd)
    pwd_command = "pwd"
    # Run the command and capture the output
    result = subprocess.run(pwd_command, shell=True, stdout=subprocess.PIPE, text=True)
    # Extract the output into a variable
    current_working_directory = result.stdout.strip()
    residue_code_csv = current_working_directory + "/pdb_files/residue_codes.csv"
    pdb_file_destination = current_working_directory + "/pdb_files/molecules"
    SmilesToPDB(smiles, name, pdb_file_destination, residue_code_csv)