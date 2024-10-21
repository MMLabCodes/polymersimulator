# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:13:57 2024

@author: danie
"""
from modules.sw_basic_functions import *
from collections import defaultdict
import os as os

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import nglview as nv
import MDAnalysis as mda
from MDAnalysis.lib import distances 
from MDAnalysis.analysis import rdf
import MDAnalysisData as data
from MDAnalysis.analysis.polymer import PersistenceLength

import pandas as pd

class master_anal():
    def __init__(self, manager, system_name, base_molecule_name, simulation_directory, poly_length=None):
        self.manager = manager
        self.system_name = system_name
        self.topology_file = self.manager.load_amber_filepaths(system_name)[0]
        self.simulation_directory = simulation_directory
        self.simulation_files = self.group_files()
        self.min_filepath = os.path.join(self.simulation_directory, self.simulation_files["min"][0])
        # It is important to note that passing a polymer length is only appropriate where the system contains polymers of the same length
        if poly_length is not None:
            self.poly_length = poly_length
            self.residue_codes = self.calculate_polymers_and_assign_residue_codes(self.min_filepath, self.poly_length)[2]
            self.poly_sel_dict = self.calculate_polymers_and_assign_residue_codes(self.min_filepath, self.poly_length)[1]
            self.number_of_polymers = len(self.poly_sel_dict)
            self.system_vol = self.base_poly_vol*self.number_of_polymers
        else:
            self.poly_length = None 
            self.residue_codes = self.extract_rescodes_and_resnums(self.min_filepath)[1]
            self.system_vol = None
        self.simulation_steps = list(self.simulation_files.keys())
        self.base_pdb = self.manager.load_pdb_filepath(base_molecule_name)
        self.base_poly_vol = vol_from_pdb(self.base_pdb)

        
    def group_files(self):
        grouped_files = defaultdict(list)

        sim_step_strings = ["1_atm", "temp_ramp_heat", "temp_ramp_cool", "min"]

        for file in os.listdir(self.simulation_directory):
            if file.endswith(('.txt', '.dcd', '.pdb')):
                base_name = os.path.splitext(file)[0]
                for string in sim_step_strings:
                    if string in base_name:
                        grouped_files[string].append(file)

        return(grouped_files)

    def extract_rescodes_and_resnums(self, pdb_file_path):
        largest_residue_number = None  # Variable to track the largest residue number
        unique_residue_codes = set()    # Set to hold unique residue codes

        with open(pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                # Parse only lines that start with "ATOM" or "HETATM"
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Extract the residue number (position 22-26)
                    residue_number = int(line[22:26].strip())
                    # Extract the residue code (position 17-20)
                    residue_code = line[17:20].strip()

                    # Update the largest residue number if this one is larger
                    if largest_residue_number is None or residue_number > largest_residue_number:
                        largest_residue_number = residue_number
                
                    # Add the residue code to the set for unique codes
                    unique_residue_codes.add(residue_code)

        return largest_residue_number, unique_residue_codes

    def calculate_polymers_and_assign_residue_codes(self, pdb_file_path, poly_length):
        # Find the largest residue number and unique residue codes
        largest_residue_number, unique_residue_codes = self.extract_rescodes_and_resnums(pdb_file_path)

        # Calculate the number of polymers
        num_polymers = largest_residue_number // poly_length

        # Create a dictionary to hold the polymer residue codes
        polymers_dict = {}

        # Assign residue codes based on the number of residues per polymer
        for i in range(num_polymers):
            # Calculate the start and end residue codes for this polymer
            start_code = i * poly_length + 1
            end_code = start_code + poly_length - 1
            polymers_dict[f'Polymer_{i + 1}'] = list(range(start_code, end_code + 1))

        return num_polymers, polymers_dict, unique_residue_codes   

class Universe():
    def __init__(self, master_anal, sim_key, traj_format=None):
        if traj_format is None:
            self.traj_format = ".pdb"
        else:
            if traj_format != ".pdb" and traj_format != ".dcd":
                print(f"{traj_format} is not supported")
                print("please enter '.pdb' or '.dcd' format.")
            else:
                self.traj_format = traj_format
        self.sim_key = sim_key
        self.masterclass = master_anal
        self.topology = self.masterclass.topology_file
        # True tells 'select_file' we are searching for the traj
        self.trajectory = os.path.join(self.masterclass.simulation_directory, self.select_file(True))
        self.universe = mda.Universe(self.topology, self.trajectory)
        self.output_filename = os.path.join(self.masterclass.simulation_directory, self.masterclass.system_name + f"_{self.sim_key}")
        # False tells 'select_file' we are searching for the data file
        self.data_file = os.path.join(self.masterclass.simulation_directory, self.select_file(False))
        self.data = pd.read_csv(self.data_file)

    def select_file(self, traj):
        if self.sim_key in self.masterclass.simulation_files:
            # Filter the files based on the specified extension
            if traj == True:
                matching_files = [filename for filename in self.masterclass.simulation_files[self.sim_key] if filename.endswith(self.traj_format)]
            if traj == False:
                matching_files = [filename for filename in self.masterclass.simulation_files[self.sim_key] if filename.endswith(".txt")]            
            
            if matching_files:
                return matching_files[0]  # Return the first matching file
            else:
                return f"No files with extension '{extension}' found for key '{self.sim_key}'."
        else:
            return f"Key '{self.sim_key}' not found in the dictionary."