# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:13:57 2024

@author: danie
"""
from modules.sw_basic_functions import *
#from modules.sw_custom_decorators import *
from collections import defaultdict
import os as os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.metrics import r2_score
import seaborn as sns
import numpy as np
import nglview as nv
import MDAnalysis as mda
from MDAnalysis.lib import distances 
from MDAnalysis.analysis import rdf
import MDAnalysisData as data
from MDAnalysis.analysis.polymer import PersistenceLength
from MDAnalysis.analysis import polymer

import pandas as pd

import warnings
warnings.filterwarnings("ignore", message="DCDReader currently makes independent timesteps")

def initialise_analysis(manager, system_name, base_molecule_name, poly_len, sim_index=0):
    '''
    This function will require an update as it only parses one type of simulation stage at the minute
    '''
    sim_avail = manager.simulations_avail(system_name)
    masterclass = master_anal(manager, system_name, base_molecule_name, sim_avail[sim_index], poly_len)
    universe = Universe(masterclass, 'temp_ramp_cool', '.dcd')
    return(universe)

class master_anal():
    def __init__(self, manager, system_name, base_molecule_name, simulation_directory, poly_length=None):
        self.manager = manager
        self.polymer_code = base_molecule_name.split("_")[0]
        self.system_name = system_name
        self.topology_file = self.manager.load_amber_filepaths(system_name)[0]
        self.simulation_directory = simulation_directory
        self.simulation_files = self.group_files()
        self.min_filepath = os.path.join(self.simulation_directory, self.simulation_files["min"][0])
        self.base_pdb = self.manager.load_pdb_filepath(base_molecule_name)
        self.base_poly_vol = estimated_volume(self.base_pdb)
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
        self.simulation_stages = list(self.simulation_files.keys())
        

        
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
    def __init__(self, master_anal, sim_stage, traj_format=None):
        if traj_format is None:
            self.traj_format = ".pdb"
        else:
            if traj_format != ".pdb" and traj_format != ".dcd":
                print(f"{traj_format} is not supported")
                print("please enter '.pdb' or '.dcd' format.")
            else:
                self.traj_format = traj_format
        self.sim_stage = sim_stage
        self.masterclass = master_anal
        self.topology = self.masterclass.topology_file
        # True tells 'select_file' we are searching for the traj
        self.trajectory = os.path.join(self.masterclass.simulation_directory, self.select_file(True))
        self.universe = mda.Universe(self.topology, self.trajectory)
        self.output_filename = os.path.join(self.masterclass.simulation_directory, self.masterclass.system_name + f"_{self.sim_stage}")
        # False tells 'select_file' we are searching for the data file
        self.data_file = os.path.join(self.masterclass.simulation_directory, self.select_file(False))
        self.data = pd.read_csv(self.data_file)

    def select_file(self, traj):
        if self.sim_stage in self.masterclass.simulation_files:
            # Filter the files based on the specified extension
            if traj == True:
                matching_files = [filename for filename in self.masterclass.simulation_files[self.sim_stage] if filename.endswith(self.traj_format)]
            if traj == False:
                matching_files = [filename for filename in self.masterclass.simulation_files[self.sim_stage] if filename.endswith(".txt")]            
            
            if matching_files:
                return matching_files[0]  # Return the first matching file
            else:
                return f"No files with extension '{extension}' found for key '{self.sim_key}'."
        else:
            return f"Key '{self.sim_stage}' not found in the dictionary."

    def select_polymer(self, polymer_name, select_rest=False):
        # Select the first polymer (Polymer_1)
        first_polymer_indices = self.masterclass.poly_sel_dict[polymer_name]

        # Convert the list of indices into a selection string
        if select_rest == False:
            selection_string = "resid " + " ".join(map(str, first_polymer_indices))
        if select_rest == True:
            selection_string = "not resid " + " ".join(map(str, first_polymer_indices))          
        selected_atoms = self.universe.select_atoms(selection_string)

        return(selected_atoms)

    def select_backbone(self, polymer_name):
        backbone_smarts = ["[O][C][C][O]", "[O][C][C][C][O]", "[O][C][C][C][C][O]", "[O][C][C][C][C][C][O]", "[O][C][C][C][C][C][C][O]",  "[O][C][C][C][C][C][C][C][O]"]
        
        # Select the first polymer (Polymer_1)
        first_polymer_indices = self.masterclass.poly_sel_dict[polymer_name]
        selection_string = "resid " + " ".join(map(str, first_polymer_indices))
        selected_atoms = self.universe.select_atoms(selection_string)

        backbone_results = []
        for smarts in backbone_smarts:
            backbone = selected_atoms.select_atoms("smarts {}".format(smarts))
            backbone_results.append(backbone)

        for result in backbone_results:
            if result.n_atoms == 0:
                pass
            else:
                return(result)

        print("No backbone identified, consider adding a backbone pattern to this functions source code")
        return(None)

        
        
        

def select_not_polymer(polymer_name):
    # Get the indices of the specified polymer from the masterclass
    polymer_indices = universe.masterclass.poly_sel_dict[polymer_name]

    # Convert the list of indices into a selection string
    # Using `not` to exclude the specified residues
    selection_string = "not resid " + " ".join(map(str, polymer_indices))

    # Use the selection string in MDAnalysis
    selected_atoms = universe.universe.select_atoms(selection_string)

    return selected_atoms

class Analysis:

    @staticmethod
    def calc_free_volume(universe, plot=None, bins=25, from_multiple=False):
        u = universe.universe
        n_frames = len(u.trajectory)
        free_volumes = []
        temperatures = universe.data["Temperature (K)"]  # Assuming temperature data is available

        for ts in u.trajectory:
            box_vectors = ts.dimensions
            volume = box_vectors[0] * box_vectors[1] * box_vectors[2]
            free_volume = volume - universe.masterclass.system_vol
            free_volumes.append(free_volume)

        avg_free_volume = sum(free_volumes) / len(free_volumes)

        if plot is None or plot is False:
            print("The returned value is an average volume across this whole simulation stage. Do not report this value if a thermal ramping was used as values will not be representative of the system at any point.")
            return avg_free_volume
        else:
            # Binning the data using np.histogram for both free volumes and temperatures
            bin_indices = np.digitize(temperatures, bins=np.linspace(200, 700, bins))
            binned_volumes = []
            binned_temperatures = []

            for i in range(1, bins + 1):
                bin_volumes = [free_volumes[j] for j in range(len(bin_indices)) if bin_indices[j] == i]
                bin_temps = [temperatures[j] for j in range(len(bin_indices)) if bin_indices[j] == i]

                if bin_volumes:
                    binned_volumes.append(np.mean(bin_volumes))
                    binned_temperatures.append(np.mean(bin_temps))

            # Return binned data if specified
            if plot == "multiple":
                return binned_temperatures, binned_volumes

            # Plot the binned data as a scatter plot
            plt.scatter(binned_temperatures, binned_volumes, marker='o', label="Data Points")

            # Fit and plot the best-fit line
            best_degree = 1
            best_r2 = -np.inf
            best_fit_line = None

            # Test polynomial fits of degree 1 (linear) to 3
            for degree in range(1, 4):
                coeffs = np.polyfit(binned_temperatures, binned_volumes, degree)
                poly = np.poly1d(coeffs)
                fit_line = poly(binned_temperatures)

                # Calculate R² score to evaluate the fit
                r2 = r2_score(binned_volumes, fit_line)
                if r2 > best_r2:
                    best_r2 = r2
                    best_degree = degree
                    best_fit_line = fit_line

            # Plot the best-fit line
            plt.plot(binned_temperatures, best_fit_line, label=f"{universe.masterclass.system_name} Best Fit (Degree {best_degree})", color='red')

            # Label the axes
            plt.xlabel('Temperature (K)')
            plt.ylabel('Free Volume (Å³)')

            # Add a title
            plt.title('Binned Free Volume vs Temperature with Best Fits')

            # Add a legend
            plt.legend()

            # Save the graph
            savepath = os.path.join(universe.masterclass.simulation_directory, "Binned_Free_Vol_vs_Temp_BestFit.png")
            plt.savefig(savepath)

            # Show the plot
            plt.show()

            return binned_temperatures, binned_volumes

    @staticmethod
    def calc_free_volume_multiple(universe_list, graph_title=None, plot=True, bins=25):
        """
        This function takes a list of universes, calls the calc_free_volume function for each, 
        and plots the results on the same graph.

        :param universe_list: List of universes
        :param plot: Boolean to indicate if you want to plot or just return the values
        :param bins: Number of bins for data binning
        :return: A list of (temperatures, volumes) for each universe
        """
        all_binned_temperatures = []
        all_binned_volumes = []
        
        # Generate colors from a colormap based on the number of universes
        colormap = cm.get_cmap('viridis', len(universe_list))  # Choose a colormap like 'viridis'

        # Loop over each universe
        for i in range(len(universe_list)):
            # Call the existing calc_free_volume function
            binned_temperatures, binned_volumes = Analysis.calc_free_volume(universe_list[i], plot="multiple", bins=bins)
            all_binned_temperatures.append(binned_temperatures)
            all_binned_volumes.append(binned_volumes)

            # Assign a unique color from the colormap
            color = colormap(i)

            # Plot the volumes for this universe using crosses as markers
            plt.scatter(
                binned_temperatures, binned_volumes, 
                color=color, label=universe_list[i].masterclass.system_name, 
                marker='x'
            )

            # Fit and plot the best-fit line for this universe
            best_degree = 1
            best_r2 = -np.inf
            best_fit_line = None

            # Test polynomial fits of degree 1 (linear) to 3
            for degree in range(1, 4):
                coeffs = np.polyfit(binned_temperatures, binned_volumes, degree)
                poly = np.poly1d(coeffs)
                fit_line = poly(binned_temperatures)

                # Calculate R² score to evaluate the fit
                r2 = r2_score(binned_volumes, fit_line)
                if r2 > best_r2:
                    best_r2 = r2
                    best_degree = degree
                    best_fit_line = fit_line

            # Plot the best-fit line
            plt.plot(binned_temperatures, best_fit_line, color=color)

        # Label the axes
        plt.xlabel('Temperature (K)')
        plt.ylabel('Free Volume (Å³)')

        # Add a title
        if graph_title == None:
            graph_title = "Free Volume vs Temperature"
            graph_filename = "Free_volume_vs_temperature"
        plt.title(graph_title)

        # Add a legend to the right side of the plot
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # Adjust layout to ensure the legend doesn't overlap with the plot
       # plt.tight_layout(rect=[0, 0, 0.85, 1])

        # Save the graph
        if graph_title != None:
            graph_filename = graph_title.replace(" ", "_")
        savepath = os.path.join(universe_list[0].masterclass.simulation_directory, graph_filename)
        plt.savefig(savepath)

        # Show plot if requested
        if plot:
            plt.show()

        # Return all data for potential further analysis
        return all_binned_temperatures, all_binned_volumes

    @staticmethod
    def plot_ROG(universe_object, atom_group, graph_filename=None, graph_title=None):
        if graph_filename == None:
            graph_filename = "_ROG_graph"
        if graph_title == None:
            graph_title = "ROG_graph"
        
        # universe = mdanalysis univers
        # atom_group = mdanalysis atom group
        # atom_group_name = str of atom_group_identifier (i.e. "_polymer", "_polymer_backbone", ect..)

        rog = []

        for ts in universe_object.universe.trajectory:
            rog.append(atom_group.radius_of_gyration())

        plt.plot(rog, linewidth=0.2)
        plt.xlabel('Time (fs)')
        plt.ylabel(r"R$_{g}$ ($\AA$)")
        plt.ylim(bottom=0)
        plt.ylim(top=20)
        plt.title(graph_title)
        if hasattr(universe_object, "output_filename"):
            # Default path is just used generally and not part of another analysis
            graph_filepath = universe_object.output_filename + graph_filename
        else:
            # Default path is just used generally and not part of another analysis
            graph_filepath = os.path.join(os.getcwd(), "ROG_graph")
        
        plt.savefig(graph_filepath)
        avg_rog = sum(rog) / len(rog)
        return(avg_rog, rog)

    @staticmethod
    def plot_ROG_5_5_array(universe_object):
        graph_title = universe_object.masterclass.system_name + "_" + universe_object.sim_stage
        avg_rogs = []
        for i in range(25):
            selected_atoms = universe_object.select_polymer("Polymer_" + str(i+1))
            rog_anal = Analysis.plot_ROG(universe_object, selected_atoms, "_5_5_array_ROG", graph_title)
            avg_rogs.append(rog_anal[0])

        return(None)

    @staticmethod
    def calc_end_to_end_dist(atom_group):
        start_pos = atom_group.positions[0]
        end_pos = atom_group.positions[-1]
        dist = distances.distance_array(start_pos, end_pos)[0][0]
        return(dist)

    @staticmethod
    def plot_end_to_end_dists_5_5_array(universe_object, plot=False):
        dists = []
        for i in range(len(universe_object.masterclass.poly_sel_dict)):
            poly = universe_object.select_polymer("Polymer_" + str(i+1))
            dist = Analysis.calc_end_to_end_dist(poly)
            dists.append(dist)
        if plot == False:
            return(dists)
        plt.figure(figsize=(8,6))
        plt.plot(range(1, len(dists) + 1), dists, marker='o', linestyle='-', color='b', label='End-to-end Distance')

        plt.xlabel("Polymer Index")
        plt.ylabel("End-to-End Distance (Å)")
        plt.title(universe_object.masterclass.system_name + "_end_to_end_dists")
        plt.legend()
        plt.grid(True)
        graph_filepath = os.path.join(universe_object.masterclass.simulation_directory, (universe_object.masterclass.system_name + "_5_5_array_end_to_end"))
        plt.savefig(graph_filepath)
        plt.show()
        return(dists)

    @staticmethod
    def plot_end_to_end_dists_5_5_array_statistics(list_of_universe_objects, graph_filename=None, graph_title=None):
        if graph_title is None:
            graph_title = "End-to-End Distances with Mean and Standard Deviation for Each System"
            graph_filepath = os.path.join(os.getcwd(), "End-to-End_graph")
        else:
            graph_filepath = os.path.join(list_of_universe_objects[0].masterclass.manager.systems_dir, graph_filename)

        system_distances = {}

        # Build the system_distances dictionary with distances as a list
        for i in range(len(list_of_universe_objects)):
            dists = Analysis.plot_end_to_end_dists_5_5_array(list_of_universe_objects[i])
            system_name = list_of_universe_objects[i].masterclass.system_name
            system_distances[system_name] = {
                "distances": dists  # Store distances as part of a dictionary
        }
    
        polymer_codes = [universe.masterclass.polymer_code for universe in list_of_universe_objects]

        # Prepare data for plotting
        system_names = list(system_distances.keys())  # List of system names for the x-axis
        distance_values = [system_distances[sys]["distances"] for sys in system_names]  # Extract distances

        # Calculate mean and standard deviation for each system
        means = [np.mean(dists) for dists in distance_values]
        std_devs = [np.std(dists) for dists in distance_values]

        # Add mean and standard deviation to the dictionary
        for i, system_name in enumerate(system_names):
            system_distances[system_name]["mean"] = means[i]
            system_distances[system_name]["std_dev"] = std_devs[i]

        # Set up the plot
        plt.figure(figsize=(10, 6))

        # Scatter plot of individual distances for each system, without adding to the legend
        for i, distances in enumerate(distance_values):
            x_vals = [i + 1] * len(distances)  # 1-based indexing for x values
            plt.scatter(x_vals, distances, alpha=0.6, label="_nolegend_")

        # Plot the means with standard deviation error bars and include in legend
        plt.errorbar(range(1, len(system_names) + 1), means, yerr=std_devs, fmt='o', color='red', 
                 capsize=5, label='Mean ± Std Dev')

        # Customize the plot
        plt.xticks(range(1, len(system_names) + 1), polymer_codes, rotation=90)
        plt.xlabel("Polymer")
        plt.ylabel("End-to-End Distance (Å)")
        plt.title(graph_title)

        # Place the legend inside the plot
        plt.legend(loc='best', borderaxespad=0.)

        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()  # Adjust layout to prevent clipping of tick labels
        plt.savefig(graph_filepath)

        # Display the plot
        plt.show()

        # Convert the system_distances dictionary to a DataFrame
        data = {
            "System": [],
            "Polymer_Code": [],
            "Mean_Distance": [],
            "Std_Dev_Distance": [],
            "All_Distances": [],
        }
        for i, system_name in enumerate(system_names):
            data["System"].append(system_name)
            data["Polymer_Code"].append(polymer_codes[i])
            data["Mean_Distance"].append(system_distances[system_name]["mean"])
            data["Std_Dev_Distance"].append(system_distances[system_name]["std_dev"])
            data["All_Distances"].append(system_distances[system_name]["distances"])

        df = pd.DataFrame(data)

        return df
        
    @staticmethod   
    def plot_PL(universe_object, atom_group, graph_filename=None, graph_title=None, plot=False):
        if graph_filename == None:
            graph_filename = "_PL_graph"
        if graph_title == None:
            graph_title = "PL_graph"
            
        pl = polymer.PersistenceLength([atom_group])
        pl.run()

        if plot == True:
            pl.plot()
            plt.title(graph_title)
            if hasattr(universe_object, "output_filename"):
                # Default path is just used generally and not part of another analysis
                graph_filepath = universe_object.output_filename + graph_filename
            else:
                # Default path is just used generally and not part of another analysis
                graph_filepath = os.path.join(os.getcwd(), "pl_graph")
        
            plt.savefig(graph_filepath)
        
        return(pl.lp, pl.lb)    
        
































