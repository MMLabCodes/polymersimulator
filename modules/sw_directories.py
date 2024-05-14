# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:18:48 2024

@author: danie
"""
import os 
import csv
'''
Polymer Simulator Directories Module

This module defines the PolymerSimulatorDirs class, which is used for initializing
and organizing directories for polymer simulation setup. It creates specific
directories such as 'python_scripts', 'pdb_files', 'molecules', 'systems', and manages
files like 'residue_codes.csv'.

Attributes:
    main_dir (str): The main directory for polymer simulation setup.

Example:
    >>> from PolymerSimulatorDirs import PolymerSimulatorDirs
    >>> polymer_dirs = PolymerSimulatorDirs('/path/to/main/dir/')
    >>> print(polymer_dirs.pdb_file_dir)
    '/path/to/main/dir/pdb_files'

Note:
    The main_dir must be provided as a valid path and must already exist.
'''
class PolymerSimulatorDirs:
    # Class for intialsing and building files and folders for polymer simulation set up and the like
    def __init__(self, main_dir): 
        """
        Initialize PolymerSimulatorDirs object.

        Args:
            main_dir (str): The main directory for polymer simulation setup.
                            Example: '/path/to/main/dir/'

        Raises:
            FileNotFoundError: If main_dir does not exist.
        """
        if not os.path.exists(main_dir):
            raise FileNotFoundError(f"The specified main directory '{main_dir}' does not exist.")       
        self.main_dir = main_dir 
        
        # Get path to python script directory and check if exists
        self.python_script_dir = os.path.join(main_dir, 'python_scripts')        
        if not os.path.exists(self.python_script_dir):
            os.makedirs(self.python_script_dir)
         
        # Get path to pdb file directory and check if exists
        self.pdb_file_dir = os.path.join(main_dir, 'pdb_files')      
        if not os.path.exists(self.pdb_file_dir):
            os.makedirs(self.pdb_file_dir)     
         
        # Get path to csv to pdb file directory and check if exists
        self.csv_to_pdb_dir = os.path.join(main_dir, 'csvs_to_pdb')
        if not os.path.exists(self.csv_to_pdb_dir):
            os.makedirs(self.csv_to_pdb_dir)
         
        # Get path to residue code csv and check if exists
        self.residue_code_csv = os.path.join(self.pdb_file_dir, 'residue_codes.csv')       
        if not os.path.exists(self.residue_code_csv):
            with open(self.residue_code_csv, 'w') as file:
                pass # Don't want to do anything with it, just generate it
        else:
            pass # Pass as the file already exists
        
        # Get path to molecules directory and check if exists - this will contain the parameters and pdb for individual molecules
        self.molecules_dir = os.path.join(self.pdb_file_dir, 'molecules')     
        if not os.path.exists(self.molecules_dir):
            os.makedirs(self.molecules_dir)
        
        # Get path to systems directory and check if exists - this will contain the parameters and pdb for systems for md simulations
        self.systems_dir = os.path.join(self.pdb_file_dir, 'systems')   
        if not os.path.exists(self.systems_dir):
            os.makedirs(self.systems_dir)

    def mol2_files_avail(self, directories):
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(directories.pdb_file_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".mol2"):
                    # Construct the full path to the .pdb file
                    mol2_file_path = os.path.join(root, file)
                    # Extract molecule name
                    mol2_file = mol2_file_path.split("/")[-1]
                    print(mol2_file)

    def load_mol2_filepath(self, directories, molecule_name):
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(directories.pdb_file_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".mol2"):
                    if molecule_name in file:
                        mol2_file_path = os.path.join(root, file)
                        return mol2_file_path 
    
    def pdb_files_avail(self):
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".pdb"):
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    print(pdb_file)

    def ac_files_avail(self):
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".ac"):
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    print(pdb_file)
    
    def load_pdb_filepath(self, molecule_name):
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                #if file.endswith(".pdb"):
                if file == (molecule_name + ".pdb"):
                    #if molecule_name in file:
                    pdb_file_path = os.path.join(root, file)
                    return pdb_file_path             
    
    def parametrized_mols_avail(self, directories):
        a = False
        for root, dirs, files in os.walk(directories.molecules_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".mol2"):
                    a = True
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    print(pdb_file)
        if a == False:
            print("No parametrized molecules.")
            
    def systems_avail(self):
        a = False
        for root, dirs, files in os.walk(self.systems_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".prmtop"):
                    a = True
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    print(pdb_file)
                if file.endswith(".rst7"):
                    a = True
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    print(pdb_file)
        if a == True:
            print("")
            print("Remember you need both .prmtop and .rst7 files to run a simulation")               
        if a == False:
            print("No parametrized molecules.")

    def retrieve_top_crds(self, system_name=None):
        if system_name == None:
            print("Please provide the name of the system you are retrieving files as follows: 'topology_file, coordinate_file = directories.retrieve_top_crds('ethane')")
            print("Change ethane for the name of the desired system")
            print("NOTE: Amber files must be generated using tleap prior to this step")
            return(None)
        for root, dirs, files in os.walk(self.systems_dir):
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                #if file.endswith(".prmtop") and system_name in file: 
                if file == (system_name + ".prmtop"):
                    # Construct the full path to the .pdb file
                    prmtop_file_path = os.path.join(root, file)
                #if file.endswith(".rst7") and system_name in file:
                if file == (system_name + ".rst7"):
                    # Construct the full path to the .pdb file
                    coord_file_path = os.path.join(root, file)   
        return(prmtop_file_path, coord_file_path)
    
    def bash_submission(self):
        pass
        # Need an if instance for whether its amber or ani
    
    @staticmethod  # Dont take instance or class - can run this without an instance   
    def unpack_csv(csv_with_mol_info):
        names = []
        smiles = []

        with open(csv_with_mol_info, 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                if len(row) >= 2:
                    names.append(row[0])
                    smiles.append(row[1])

        return names, smiles
        
    def retrieve_polymeric_rescodes(self, directories, molecule_name):
        with open(directories.residue_code_csv, 'r') as file:
            head_code, mainchain_code, tail_code = None, None, None
            for line in file:
                parts = (line.strip().split('\t'))[0]
                if len(parts.split(",")) >= 3 and parts.split(",")[0] == ("head_" + molecule_name):
                    head_code = parts.split(",")[2]
                if len(parts.split(",")) >= 3 and parts.split(",")[0] == ("mainchain_" + molecule_name):
                    mainchain_code = parts.split(",")[2]
                if len(parts.split(",")) >= 3 and parts.split(",")[0] == ("tail_" + molecule_name):
                    tail_code = parts.split(",")[2]    
        return(head_code, mainchain_code, tail_code)
    # Need to add stuff for md results and the like
    # Need to add a way to import the bash submission scripts - will get to this at some point
    