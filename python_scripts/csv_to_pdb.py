# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 16:51:16 2023

@author: danie
"""
import csv
import sys as sys
import os
import subprocess
from sw_openmm import SmilesToPDB, SmilesToPDB_GenerateRescode
from sw_directories import PolymerSimulatorDirs
'''
Smiles to PDB Conversion Script - using a csv file

This script reads molecule names and SMILES strings from a CSV file and generates
PDB files with optional residue codes for each molecule using the sw_openmm module.
It utilizes the PolymerSimulatorDirs class from the sw_directories module to manage
directory paths for generation of pdb files from SMILES strings.

Usage:
    python3 script_name.py <csv_file>

    <csv_file>: CSV filename with names and SMILES strings of each molecule.

Example:
    python3 script_name.py molecules.csv

Note:
    - Ensure that the required modules are installed before running the script.
    - The CSV file should have two columns: 'Name' and 'SMILES'.
    - The generated PDB files will be saved in the 'pdb_files' directory.
'''
if __name__ == "__main__":
    
    def print_error_message():
        """
        Print an error message for incorrect command-line arguments.
        """
        print("Error: Incorrect number of arguments.")
        print("Usage: python3 script_name.py <csv_filename>")
        print("Example: python3 script_name.py alkanes.csv")
        print("  <csv_file>: CSV file with names and smilestrings of each molecule.")
        sys.exit(1)
        
    # Check the number of arguments
    if len(sys.argv) != 2:
        print_error_message()

    # Bash command to print the current working directory (pwd)
    pwd_command = "pwd"

    # Run the command and capture the output
    result = subprocess.run(pwd_command, shell=True, stdout=subprocess.PIPE, text=True)

    # Extract the output into a variable
    current_working_directory = result.stdout.strip()

    # Now get our directory paths from the current_working_directory-
    directories = PolymerSimulatorDirs(current_working_directory)

    # Unpack smiles and names from the csv file
    csv_filename = os.path.join(directories.csv_to_pdb_dir, sys.argv[1])
    names_list, smiles_list = PolymerSimulatorDirs.unpack_csv(csv_filename)
    
    # Now the names_list and smiles_list contain the respective lists from the CSV file.
    for i in range(len(smiles_list)):
        SmilesToPDB_GenerateRescode(smiles_list[i], names_list[i], directories.pdb_file_dir, directories.residue_code_csv)
        



