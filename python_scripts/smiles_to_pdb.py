# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:51:58 2024

@author: danie
"""
import sys
import subprocess
from sw_openmm import SmilesToPDB, SmilesToPDB_GenerateRescode
from sw_directories import PolymerSimulatorDirs
'''
"""
Smiles to PDB Conversion Script

This script converts a SMILES string to a PDB file using the sw_openmm library.
It provides a command-line interface for specifying the name, SMILES string, and
whether to generate residue codes.

Usage:
    python3 script.py <name> <smiles> <gen_code>

Example:
    python3 smiles_to_pdb.py 'methane' 'C' 'True'

Parameters:
    - name (str): The name for the molecule.
    - smiles (str): The SMILES string representing the molecule.
    - gen_code (str): A string indicating whether to generate residue codes.
                     Accepted values: 'True' or 'False'.

Note:
    Ensure that the required Python libraries, sw_openmm and sw_directories,
    are installed before running the script.
'''
if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Error: Insufficient command-line arguments.")
        print("Usage: python3 script.py <name> <smiles> <gen_code>")
        print("Example: python3 smiles_to_pdb.py 'methane' 'C' 'True'")
        sys.exit(1)

    name = sys.argv[1]
    smiles = sys.argv[2]
    gen_code = sys.argv[3]
    
    # Bash command to print the current working directory (pwd)
    pwd_command = "pwd"
    
    # Run the command and capture the output
    result = subprocess.run(pwd_command, shell=True, stdout=subprocess.PIPE, text=True)
    
    # Extract the output into a variable
    current_working_directory = result.stdout.strip()
    
    # Now get our directory paths from the current_working_directory-
    directories = PolymerSimulatorDirs(current_working_directory)
    
    # Perform Smiles to PDB conversion based on gen_code value
    if gen_code == 'True':
        SmilesToPDB_GenerateRescode(smiles, name, directories.pdb_file_dir, directories.residue_code_csv)
    if gen_code == 'False':
        SmilesToPDB(smiles, directories.pdb_file_dir)