# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:13:57 2024

@author: danie
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import nglview as nv
import MDAnalysis as mda
from MDAnalysis.lib import distances 
from MDAnalysis.analysis import rdf
import MDAnalysisData as data
import re
import os

class Analysis():
    cached_ROG_data = None
    cached_ROG_average = None
    cached_COG_data = None
    cached_COG_average = None
    
    
    def __init__(self):
        pass

    @staticmethod
    def plot_ROG(atom_group, atom_group_name, graph_filepath=None):
        # NOTE: there is an issue when trying to save the graph as a png

        rog = [] # List to storeresults
        #Os = u.select_atoms('type os') # Select atom types

        # Iterate over each frame and collect ROG
        for ts in self.u.trajectory:
            rog.append(atom_group.radius_of_gyration())

        # Plot the frequencies of the residue names
        plt.plot(rog)
        plt.figure(figsize=(10, 6))
        plt.xlabel('frame')
        plt.ylabel(r"R$_{g}$ ($\AA$)")
       # plt.title()
        # Save the plot to a file
        if graph_filepath == None:
            # Default path is just used generally and not part of another analysis
            graph_filepath = os.path.join(os.getcwd(), "ROG_graph")
            plt.savefig(graph_filepath)
        else:
            plt.savefig(graph_filepath)
        plt.close()

        Analysis.cached_ROG_data = rog
        Analysis.cached_ROG_average = sum(rog) / len(rog)
        return(None)

    @staticmethod
    def plot_COG(atom_group, atom_group_name, graph_filepath=None):
        cog = [] # List to storeresults
        #Os = u.select_atoms('type os') # Select atom types

        # Iterate over each frame and collect ROG
        for ts in self.u.trajectory:
            cog.append(atom_group.center_of_mass())

        cog = np.array(cog)
        # Plot the frequencies of the residue names
        # Now let's plot each coordinate separately
        plt.plot(cog[:, 0], label='x')
        plt.plot(cog[:, 1], label='y')
        plt.plot(cog[:, 2], label='z')
        plt.figure(figsize=(10, 6))
        plt.xlabel('frame')
        plt.ylabel(r"R$_{g}$ ($\AA$)")
       # plt.title()
        # Save the plot to a file
        if graph_filepath == None:
            # Default path is just used generally and not part of another analysis
            graph_filepath = os.path.join(os.getcwd(), "ROG_graph")
            plt.savefig(graph_filepath)
        else:
            plt.savefig(graph_filepath)
        plt.close()

        Analysis.cached_COG_data = cog
        Analysis.cached_COG_average = (sum(cog[:, 0]) / len(cog[:, 0])), (sum(cog[:, 1]) / len(cog[:, 1])), (sum(cog[:, 2]) / len(cog[:, 2]))
        return(None)

class SinglePolyAnalysis(Analysis):
    # Class for intialsing and building files and folders for polymer simulation set up and the like
    def __init__(self, topology_file, trajectory_file):
        print("This function is for single solvated polymers")
        print("Do you wish to continue?")
        print("Enter: y/n")
        print("")
        confirmation = input()
        if confirmation == "n":
            print("")
            print("Instance creation aborted by user.")
            return      
        if confirmation == "y":
            self.u = mda.Universe(topology_file, trajectory_file)
            self.topology_file = topology_file
            self.trajectory_file = trajectory_file
            self.output_filepath = os.path.dirname(self.trajectory_file)

    def extract_filename(self, filepath):
        # Extract the filename from the full path
        filename_with_extension = filepath.split('/')[-1]
    
        # Remove the date and time (pattern: YYYY-MM-DD_HHMMSS)
        filename_without_date = re.sub(r'_\d{4}-\d{2}-\d{2}_\d{6}', '', filename_with_extension)
    
        # Remove the file extension
        filename_without_extension = filename_without_date.split('.')[0]
    
        return(filename_without_extension)
        
    def select_polymers(self):
        unique_residues = []
        for residue in self.u.residues:
            if residue.resname not in unique_residues:
                unique_residues.append(residue.resname)
                polymer_atom_group = self.u.select_atoms('not resname WAT')
    
        return(polymer_atom_group)

    def plot_ROG(self, atom_group, atom_group_name):
        graph_filepath = os.path.join(self.output_filepath, (self.extract_filename(self.trajectory_file) + "_" + atom_group_name + "_ROG_graph"))
        super().plot_ROG(self, atom_group, atom_group_name, graph_filepath)

