# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:13:57 2024

@author: danie
"""
from modules.sw_basic_functions import *
from modules.sw_complex_fluid_models import *
#from modules.sw_custom_decorators import *
from collections import defaultdict
import os as os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
from scipy.integrate import quad
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
import pwlf

import warnings
warnings.filterwarnings("ignore", message="DCDReader currently makes independent timesteps")

def initialise_poly_analysis(manager, system_name, base_molecule_name, poly_len=10, sim_stage_name=None, sim_index=0):
    """
    Initializes polymer analysis by setting up a molecular dynamics universe.
    
    Parameters:
    manager : object
        Manages simulation files and data.
    system_name : str
        Name of the system being analyzed.
    base_molecule_name : str
        Name of the base polymer molecule.
    poly_len : int
        Length of the polymer chains.
    sim_index : int, optional
        Index of the simulation stage to use (default is 0).
    
    Returns:
    poly_Universe
        An initialized polymer universe for analysis.
    """
    if sim_stage_name == None:
        print("Please provide a simulation stage name - will be output filename tag assigned to specific stage of simulation.")
        return()
    sim_avail = manager.simulations_avail(system_name)
    masterclass = master_poly_anal(manager, system_name, base_molecule_name, sim_avail[sim_index], poly_len)
    return poly_Universe(masterclass, sim_stage_name, '.dcd')

def initialise_bio_oil_analysis(manager, system_name, sim_index=0, sim_step=None):
    """
    Initializes bio-oil analysis by setting up a molecular dynamics universe.
    
    Parameters:
    manager : object
        Manages simulation files and data.
    system_name : str
        Name of the system being analyzed.
    sim_index : int, optional
        Index of the simulation stage to use (default is 0).
    
    Returns:
    bio_oil_Universe
        An initialized bio-oil universe for analysis.
    """
    sim_avail = manager.simulations_avail(system_name)
    masterclass = master_bio_oil_anal(manager, system_name, sim_avail[sim_index])
    residue_codes = complex_fluid_model_builder.extract_unique_rescodes(masterclass.min_filepath)
    molecule_names = complex_fluid_model_builder.find_matching_molecules(residue_codes, complex_fluid_model_builder.load_molecule_list(manager))
    molecule_dictionary = dict(zip(molecule_names, residue_codes))
    if sim_step == None:
        sim_step = 'prod'
    return bio_oil_Universe(masterclass, sim_step, molecule_dictionary)

class initialise:
    """Base class for initializing and grouping simulation files."""
    
    def __init__(self):
        pass

    def group_files(self):
        """
        Groups simulation files based on predefined simulation steps.
        
        Returns:
        dict
            Dictionary mapping simulation steps to lists of corresponding files.
        """
        grouped_files = defaultdict(list)
        sim_step_strings = ["1_atm", "temp_ramp_heat", "temp_ramp_cool", "min", "prod", "cooling_NPT_cool"]

        for file in os.listdir(self.simulation_directory):
            if file.endswith(('.txt', '.dcd', '.pdb')):
                base_name = os.path.splitext(file)[0]
                for string in sim_step_strings:
                    if string in base_name:
                        grouped_files[string].append(file)
        return grouped_files

class master_bio_oil_anal(initialise):
    """
    Class for managing bio-oil molecular dynamics simulations.

    Attributes:
        manager: Object responsible for handling file operations.
        system_name (str): Name of the system being analyzed.
        topology_file (str): Path to the topology file.
        simulation_directory (str): Directory where simulation files are stored.
        simulation_files (dict): Grouped simulation files categorized by stage.
        min_filepath (str): Filepath for the minimization stage.
        simulation_stages (list): List of simulation stages.
    """

    def __init__(self, manager, system_name, simulation_directory):
        self.manager = manager
        self.system_name = system_name
        self.topology_file = self.manager.load_amber_filepaths(system_name)[0]
        self.simulation_directory = simulation_directory
        self.simulation_files = self.group_files()
        self.min_filepath = os.path.join(self.simulation_directory, self.simulation_files["min"][0])
        self.simulation_stages = list(self.simulation_files.keys())

class master_poly_anal(initialise):
    """
    Class for analyzing polymer molecular dynamics simulations.

    Attributes:
        manager: Object responsible for handling file operations.
        system_name (str): Name of the system being analyzed.
        base_molecule_name (str): Base molecule identifier.
        topology_file (str): Path to the topology file.
        simulation_directory (str): Directory where simulation files are stored.
        extracted_poly_coords_dir (str): Directory for extracted polymer coordinates.
        poly_dft_input_dir (str): Directory for polymer DFT input files.
        poly_dft_output_dir (str): Directory for polymer DFT output files.
        simulation_files (dict): Grouped simulation files categorized by stage.
        min_filepath (str): Filepath for the minimization stage.
        base_pdb (str): Path to the base polymer structure file.
        base_poly_vol (float): Estimated volume of the base polymer.
        poly_length (int, optional): Length of polymer chains (if applicable).
        residue_codes (set): Unique residue codes in the system.
        system_vol (float, optional): Total system volume if polymer length is provided.
        simulation_stages (list): List of simulation stages.
    """

    def __init__(self, manager, system_name, base_molecule_name, simulation_directory, poly_length=None):
        """
        Initializes the master_poly_anal class.

        Args:
            manager: Object responsible for handling file operations.
            system_name (str): Name of the system being analyzed.
            base_molecule_name (str): Base molecule identifier.
            simulation_directory (str): Path to the simulation directory.
            poly_length (int, optional): Length of polymer chains (default: None).
        """
        self.manager = manager
        self.polymer_code = base_molecule_name.split("_")[0]
        self.system_name = system_name
        self.base_molecule_name = base_molecule_name
        self.topology_file = self.manager.load_amber_filepaths(system_name)[0]
        self.simulation_directory = simulation_directory

        # Create necessary directories if they don't exist
        for dir_name in ["extracted_poly_coords", "poly_DFT_inputs", "poly_DFT_outputs"]:
            dir_path = os.path.join(self.simulation_directory, dir_name)
            os.makedirs(dir_path, exist_ok=True)

        self.simulation_files = self.group_files()
        self.min_filepath = os.path.join(self.simulation_directory, self.simulation_files["min"][0])
        self.base_pdb = self.manager.load_pdb_filepath(base_molecule_name)
        self.base_poly_vol = estimated_volume(self.base_pdb)

        # Handle polymer length-dependent attributes
        if poly_length is not None:
            self.poly_length = poly_length
            num_polymers, self.poly_sel_dict, self.residue_codes = self.calculate_polymers_and_assign_residue_codes(
                self.min_filepath, self.poly_length
            )
            self.number_of_polymers = len(self.poly_sel_dict)
            self.system_vol = self.base_poly_vol * self.number_of_polymers
        else:
            self.poly_length = None
            _, self.residue_codes = self.extract_rescodes_and_resnums(self.min_filepath)
            self.system_vol = None

        self.simulation_stages = list(self.simulation_files.keys())

    def extract_rescodes_and_resnums(self, pdb_file_path):
        """
        Extracts unique residue codes and finds the highest residue number from a PDB file.

        Args:
            pdb_file_path (str): Path to the PDB file.

        Returns:
            tuple: (largest_residue_number, unique_residue_codes)
                - largest_residue_number (int): The highest residue number found.
                - unique_residue_codes (set): Set of unique residue codes.
        """
        largest_residue_number = None
        unique_residue_codes = set()

        with open(pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith(("ATOM", "HETATM")):
                    residue_number = int(line[22:26].strip())  # Extract residue number
                    residue_code = line[17:20].strip()  # Extract residue code

                    # Update the largest residue number
                    if largest_residue_number is None or residue_number > largest_residue_number:
                        largest_residue_number = residue_number
                
                    unique_residue_codes.add(residue_code)

        return largest_residue_number, unique_residue_codes

    def calculate_polymers_and_assign_residue_codes(self, pdb_file_path, poly_length):
        """
        Assigns residue codes to polymers and determines polymer count.

        Args:
            pdb_file_path (str): Path to the PDB file.
            poly_length (int): Length of polymer chains.

        Returns:
            tuple: (num_polymers, polymers_dict, unique_residue_codes)
                - num_polymers (int): Number of polymers in the system.
                - polymers_dict (dict): Dictionary mapping polymer names to residue indices.
                - unique_residue_codes (set): Set of unique residue codes.
        """
        largest_residue_number, unique_residue_codes = self.extract_rescodes_and_resnums(pdb_file_path)
        num_polymers = largest_residue_number // poly_length

        # Assign residue indices to polymers
        polymers_dict = {
            f'Polymer_{i + 1}': list(range(i * poly_length + 1, (i + 1) * poly_length + 1))
            for i in range(num_polymers)
        }

        return num_polymers, polymers_dict, unique_residue_codes   

class Universe:
    """
    A class to represent a molecular dynamics universe, handling trajectories, topologies, and data files.

    Attributes:
        masterclass: An instance of the master analysis class containing simulation details.
        sim_stage (str): The stage of the simulation (e.g., 'prod', 'min', etc.).
        traj_format (str): The format of the trajectory file ('.pdb' or '.dcd').
        topology (str): Path to the topology file.
        trajectory (str): Path to the trajectory file.
        universe (MDAnalysis.Universe): MDAnalysis Universe object containing topology and trajectory data.
        output_filename (str): Output file path for storing results.
        data_file (str): Path to the associated data file.
        data (pd.DataFrame): Loaded simulation data.
    """

    def __init__(self, master_anal, sim_stage, traj_format=None):
        """
        Initializes the Universe object.

        Args:
            master_anal: An instance of the master analysis class.
            sim_stage (str): The stage of the simulation.
            traj_format (str, optional): The trajectory file format (default is '.pdb').
        """
        self.traj_format = traj_format if traj_format in [".pdb", ".dcd"] else ".pdb"

        if traj_format not in [None, ".pdb", ".dcd"]:
            print(f"Warning: {traj_format} is not supported. Using default format '.pdb'.")

        self.sim_stage = sim_stage
        self.masterclass = master_anal
        self.topology = self.masterclass.topology_file

        # Load trajectory and data files
        self.trajectory = os.path.join(self.masterclass.simulation_directory, self.select_file(traj=True))
        self.universe = mda.Universe(self.topology, self.trajectory)
        self.output_filename = os.path.join(self.masterclass.simulation_directory, f"{self.masterclass.system_name}_{self.sim_stage}")

        self.data_file = os.path.join(self.masterclass.simulation_directory, self.select_file(traj=False))
        self.data = pd.read_csv(self.data_file)

    def select_file(self, traj):
        """
        Selects the appropriate trajectory or data file based on the simulation stage.

        Args:
            traj (bool): If True, selects a trajectory file; if False, selects a data file.

        Returns:
            str: The selected filename if found, otherwise an error message.
        """
        if self.sim_stage not in self.masterclass.simulation_files:
            return f"Error: Key '{self.sim_stage}' not found in the simulation files dictionary."

        # Determine the file extension based on the file type
        ext = self.traj_format if traj else ".txt"
        matching_files = [f for f in self.masterclass.simulation_files[self.sim_stage] if f.endswith(ext)]

        return matching_files[0] if matching_files else f"Error: No files with extension '{ext}' found for key '{self.sim_stage}'."

class poly_Universe(Universe):
    """
    A subclass of Universe specifically designed for handling polymer simulations.

    Attributes:
        masterclass: An instance of the master analysis class containing polymer simulation details.
        sim_stage (str): The stage of the simulation (e.g., 'prod', 'min', etc.).
        traj_format (str): The format of the trajectory file ('.pdb' or '.dcd').
        topology (str): Path to the topology file.
        trajectory (str): Path to the trajectory file.
        universe (MDAnalysis.Universe): MDAnalysis Universe object containing topology and trajectory data.
        output_filename (str): Output file path for storing results.
        data_file (str): Path to the associated data file.
        data (pd.DataFrame): Loaded simulation data.
    """

    def __init__(self, master_anal, sim_stage, traj_format=None):
        """
        Initializes the poly_Universe object.

        Args:
            master_anal: An instance of the master analysis class.
            sim_stage (str): The stage of the simulation.
            traj_format (str, optional): The trajectory file format (default is '.pdb').
        """
        self.traj_format = traj_format if traj_format in [".pdb", ".dcd"] else ".pdb"

        if traj_format not in [None, ".pdb", ".dcd"]:
            print(f"Warning: {traj_format} is not supported. Using default format '.pdb'.")

        self.sim_stage = sim_stage
        self.masterclass = master_anal
        self.topology = self.masterclass.topology_file

        # Load trajectory and data files
        self.trajectory = os.path.join(self.masterclass.simulation_directory, self.select_file(traj=True))
        self.universe = mda.Universe(self.topology, self.trajectory)
        self.output_filename = os.path.join(self.masterclass.simulation_directory, f"{self.masterclass.system_name}_{self.sim_stage}")

        self.data_file = os.path.join(self.masterclass.simulation_directory, self.select_file(traj=False))
        self.data = pd.read_csv(self.data_file)

    def select_polymer(self, polymer_name, select_rest=False):
        """
        Selects atoms corresponding to a specified polymer.

        Args:
            polymer_name (str): The name of the polymer to select.
            select_rest (bool, optional): If True, selects all atoms *except* the specified polymer (default: False).

        Returns:
            MDAnalysis.AtomGroup: Selected atoms based on polymer selection criteria.
        """
        first_polymer_indices = self.masterclass.poly_sel_dict[polymer_name]

        # Construct the selection string
        selection_string = ("not " if select_rest else "") + "resid " + " ".join(map(str, first_polymer_indices))
        return self.universe.select_atoms(selection_string)

    def select_backbone(self, polymer_name):
        """
        Selects the polymer backbone based on predefined SMARTS patterns.

        Args:
            polymer_name (str): The name of the polymer to select.

        Returns:
            MDAnalysis.AtomGroup or None: The identified backbone atoms, or None if no match is found.
        """
        backbone_smarts = [
            "[O][C][C][O]", "[O][C][C][C][O]", "[O][C][C][C][C][O]",
            "[O][C][C][C][C][C][O]", "[O][C][C][C][C][C][C][O]", "[O][C][C][C][C][C][C][C][O]"
        ]

        first_polymer_indices = self.masterclass.poly_sel_dict[polymer_name]
        selection_string = "resid " + " ".join(map(str, first_polymer_indices))
        selected_atoms = self.universe.select_atoms(selection_string)

        # Attempt to match the backbone using predefined SMARTS patterns
        for smarts in backbone_smarts:
            backbone = selected_atoms.select_atoms(f"smarts {smarts}")
            if backbone.n_atoms > 0:
                return backbone

        print("No backbone identified. Consider adding a new SMARTS pattern to this function.")
        return None

class bio_oil_Universe(Universe):
    """
    A subclass of Universe specifically designed for handling bio-oil molecular dynamics simulations.

    Attributes:
        masterclass: An instance of the master analysis class containing bio-oil simulation details.
        sim_stage (str): The stage of the simulation (e.g., 'prod', 'min', etc.).
        traj_format (str): The format of the trajectory file ('.pdb' or '.dcd').
        topology (str): Path to the topology file.
        trajectory (str): Path to the trajectory file.
        universe (MDAnalysis.Universe): MDAnalysis Universe object containing topology and trajectory data.
        output_filename (str): Output file path for storing results.
        data_file (str): Path to the associated data file.
        data (pd.DataFrame): Loaded simulation data.
        molecule_dictionary (dict): Dictionary mapping molecule names to their residue codes.
    """

    def __init__(self, master_anal, sim_stage, molecule_dictionary, traj_format=None):
        """
        Initializes the bio_oil_Universe object.

        Args:
            master_anal: An instance of the master analysis class.
            sim_stage (str): The stage of the simulation.
            molecule_dictionary (dict): Dictionary of molecule names and corresponding residue codes.
            traj_format (str, optional): The trajectory file format (default is '.dcd').
        """
        self.traj_format = traj_format if traj_format in [".pdb", ".dcd"] else ".dcd"

        if traj_format not in [None, ".pdb", ".dcd"]:
            print(f"Warning: {traj_format} is not supported. Using default format '.dcd'.")

        self.sim_stage = sim_stage
        self.masterclass = master_anal
        self.topology = self.masterclass.topology_file

        # Load trajectory and data files
        self.trajectory = os.path.join(self.masterclass.simulation_directory, self.select_file(traj=True))
        self.universe = mda.Universe(self.topology, self.trajectory)
        self.output_filename = os.path.join(self.masterclass.simulation_directory, f"{self.masterclass.system_name}_{self.sim_stage}")

        self.data_file = os.path.join(self.masterclass.simulation_directory, self.select_file(traj=False))
        self.data = pd.read_csv(self.data_file)

        # Store molecule dictionary for reference
        self.molecule_dictionary = molecule_dictionary
        
        

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
        """
        Calculate the free volume of a system over a trajectory and optionally analyze its dependence on temperature.

        Parameters:
        -----------
        universe : MDAnalysis.Universe
            The MDAnalysis Universe object containing trajectory data and system parameters.
        plot : {None, False, "multiple"}, optional
            - If None or False (default), returns the average free volume over the trajectory.
            - If "multiple", returns binned temperature and free volume data without plotting.
            - Otherwise, generates and saves a scatter plot of binned free volume vs. temperature.
        bins : int, optional
            The number of bins for temperature-based data grouping (default: 25).
        from_multiple : bool, optional
            Reserved parameter, currently unused.

        Returns:
        --------
        float
            If `plot` is None or False, returns the average free volume over the trajectory.
        tuple of (list, list)
            If `plot` is "multiple" or plotting is enabled, returns two lists:
            - `binned_temperatures`: The average temperature in each bin.
            - `binned_volumes`: The corresponding average free volume in each bin.

        Notes:
        ------
        - The free volume is computed as the difference between the total simulation box volume and 
          the system's intrinsic volume (`universe.masterclass.system_vol`).
        - If a thermal ramping protocol was used, reporting the average free volume may be misleading.
        - When plotting is enabled, the function fits polynomial curves (degree 1 to 3) to the binned data 
          and selects the best fit based on the R² score.
        - The scatter plot and best-fit line are saved as `"Binned_Free_Vol_vs_Temp_BestFit.png"` 
          in the simulation directory.

        Example Usage:
        --------------
        >>> avg_volume = calc_free_volume(universe)
        >>> temps, volumes = calc_free_volume(universe, plot="multiple")
        >>> calc_free_volume(universe, plot=True)  # Generates and saves a plot

        """
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

        return(avg_rogs)

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

    @staticmethod
    def calculate_PL_5_5_array(universe_object):
        pls = []
        pbs = []
        for i in range(len(universe_object.masterclass.poly_sel_dict)):
            backbone = universe_object.select_backbone("Polymer_" + str(i+1))
            lp, lb = Analysis.plot_PL(universe_object, backbone)
            pls.append(lp)
            pbs.append(lb)
    
        return(pls, pbs)
    
    @staticmethod
    def plot_PL_5_5_array(universe_object, plot=False):
        pls, pbs = Analysis.calculate_PL_5_5_array(universe_object)
        if not plot:
            return pls, pbs
    
        # Plot for persistence lengths (PLs)
        plt.figure(figsize=(8, 6))
        plt.plot(range(1, len(pls) + 1), pls, marker='o', linestyle='-', color='b', label='End-to-End Distance')
        plt.xlabel("Polymer Index")
        plt.ylabel("Persistence Length (Å)")
        plt.title(universe_object.masterclass.system_name + "_persistence_lengths")
        plt.legend()
        plt.grid(True)
        graph_filepath_pls = os.path.join(
            universe_object.masterclass.simulation_directory,
            universe_object.masterclass.system_name + "_5_5_array_persistence_lengths.png"
        )
        plt.savefig(graph_filepath_pls)
        plt.show()
    
        # Plot for bond lengths (PBs)
        plt.figure(figsize=(8, 6))
        plt.plot(range(1, len(pbs) + 1), pbs, marker='o', linestyle='-', color='g', label='Bond Lengths')
        plt.xlabel("Polymer Index")
        plt.ylabel("Bond Length (Å)")
        plt.title(universe_object.masterclass.system_name + "_bond_lengths")
        plt.legend()
        plt.grid(True)
        graph_filepath_pbs = os.path.join(
            universe_object.masterclass.simulation_directory,
            universe_object.masterclass.system_name + "_5_5_array_bond_lengths.png"
        )
        plt.savefig(graph_filepath_pbs)
        plt.show()
    
        return pls, pbs

    @staticmethod
    def plot_expansion_coeff(universe_object, bin_params=None, fit=True, plot=None):
        if bin_params is None:
            print("Please specify the bin parameters for your system: [start_temp, target_temp, temp_increments]")
            return None, None

        # Set title
        graph_title = f"{universe_object.masterclass.system_name}_{universe_object.sim_stage} vol. vs temp."

        # Define dataframe
        df = universe_object.data

        # Define bin edges
        bin_width = bin_params[2]  # Temperature bin size
        bins = np.arange(bin_params[0], bin_params[1] + bin_width, bin_width)

        # Assign temperature values to bins
        df["Temp_Bin"] = pd.cut(df["Temperature (K)"], bins, labels=(bins[:-1] + bin_width / 2))

        # Group by bins and calculate average box volume
        binned_data = df.groupby("Temp_Bin")["Box Volume (nm^3)"].mean().reset_index()
        binned_data["Temp_Bin"] = binned_data["Temp_Bin"].astype(float)

        # Extract binned temperature and volume data
        T_bin = binned_data["Temp_Bin"].values
        V_bin = binned_data["Box Volume (nm^3)"].values

        fit_equations = {}
        params = None

        if fit:
            if plot == True:
                plt.figure(figsize=(8, 6))
                plt.scatter(T_bin, V_bin, marker="o", color="b", label="V(T)")
        
            # Fit model to data
            params, _ = curve_fit(volume_model, T_bin, V_bin)
            a, b, c = params

            # Generate smooth data for plotting
            T_smooth = np.linspace(min(T_bin), max(T_bin), 200)
            V_smooth = volume_model(T_smooth, *params)

            # Compute smoothed dV/dT and thermal expansion coefficient α(T)
            dVdT_smooth = b + 2*c*T_smooth
            alpha_smooth = (dVdT_smooth / V_smooth) * 1e4  # Scale factor for visualization

            # Compute α(T) at each temperature bin
            dVdT_bin = b + 2*c*T_bin
            alpha_bin = (dVdT_bin / V_bin) * 1e4

            # Store fit equations
            fit_equations["Volume Fit"] = f"V(T) = {a:.4f} + {b:.4f}T + {c:.4f}T²"
            fit_equations["Expansion Coefficient Fit"] = f"α(T) = ({b:.4f} + 2*{c:.4f}T) / V(T) * 10⁴"

            if plot == True:
                # Plot fitted volume curve
                plt.plot(T_smooth, V_smooth, linestyle="--", color="r", label="Fitted V(T)")
                plt.xlabel("Temperature (K)", fontsize=16)
                plt.ylabel("Average Box Volume (nm³)", fontsize=16)
                plt.title(graph_title, fontsize=18)
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                plt.grid(True)
                plt.legend(fontsize=14)
                plt.show()

                # Plot thermal expansion coefficient
                plt.figure(figsize=(8, 6))
                plt.scatter(T_bin, alpha_bin, color="b", label="α(T)")
                plt.plot(T_smooth, alpha_smooth, linestyle="--", color="g", label="Fitted α(T)")
                plt.xlabel("Temperature (K)", fontsize=16)
                plt.ylabel(r"Thermal Expansion Coefficient α ($\times 10^{-4}$ K⁻¹)", fontsize=16)
                plt.title("Thermal Expansion Coefficient vs Temperature", fontsize=18)
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                plt.grid(True)
                plt.legend(fontsize=14)
                plt.show()

        return fit_equations, params

    @staticmethod
    def predict_expansion(params, V1, T1, T2):
        if params is None:
            print("No fitted parameters available. Please fit the model first.")
            return None
    
        a, b, c = params
    
        # Define the fitted alpha(T) function
        def alpha_T_func(T):
            return (b + 2 * c * T) / volume_model(T, a, b, c)  # α(T) = (dV/dT) / V(T)
    
        # Integrate α(T) over the given range
        integral_alpha, _ = quad(alpha_T_func, T1, T2)
    
        # Compute total volume change ΔV
        delta_V = V1 * integral_alpha  # Integrated expansion effect
    
        return delta_V

    @staticmethod
    def plot_tg(universe, min_temp=None, max_temp=None, quench=None):
        df = universe.data
    
        # Create bins for temperature with 20 K intervals
        bins = pd.cut(df['Temperature (K)'], bins=range(min_temp, max_temp+quench, quench))
    
        # Compute the average density per bin
        average_density = df.groupby(bins)['Density (g/mL)'].mean()
        temp_midpoints = average_density.index.categories.mid
    
        # Use pwlf for piecewise linear fitting with automatic segment selection
        pwlf_model = pwlf.PiecewiseLinFit(temp_midpoints, average_density)
    
        # Determine the optimal number of segments using AIC (Akaike Information Criterion)
        max_segments = 5  # Limit to prevent overfitting
        best_segments = 1
        best_aic = float('inf')
        best_breaks = None
    
        for n in range(2, max_segments + 1):
            breaks = pwlf_model.fit(n)
            aic = pwlf_model.ssr + 2 * n  # AIC approximation
            if aic < best_aic:
                best_aic = aic
                best_segments = n
                best_breaks = breaks
    
        # Fit using the optimal number of segments
        pwlf_model.fit(best_segments)
    
        # Get the equations of the linear segments
        equations = []
        for i in range(best_segments):
            slope, intercept = pwlf_model.slopes[i], pwlf_model.intercepts[i]
            equations.append(f'y = {slope:.3f}x + {intercept:.3f}')
    
        # Print segment equations and breakpoints
        for i, eq in enumerate(equations):
            print(f'Segment {i + 1}: {eq}')
        print(f'Breakpoints at: {best_breaks[1:-1]} K')
    
        # Plot the data
        plt.plot(temp_midpoints, average_density, 'o', label='Average Density')
        predicted_density = pwlf_model.predict(temp_midpoints)
        plt.plot(temp_midpoints, predicted_density, 'r-', label=f'Piecewise Fit')
    
        # Plot the breakpoints
        for bp in best_breaks[1:-1]:  # Ignore first and last since they are dataset bounds
            plt.axvline(bp, linestyle='--', color='black', alpha=0.6)
            plt.text(bp, min(average_density), f'{bp:.2f} K', rotation=90, verticalalignment='bottom', fontsize=14)
    
        # Labels and formatting
        plt.xlabel('Temperature (K)', fontsize=14)
        plt.ylabel('Average Density (g/mL)', fontsize=14)
        plt.title(f'Density vs temperature for {universe.masterclass.system_name}', fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend()
        plt.grid(True)
    
        # Show the plot
        plt.show()
        return(float(f'{best_breaks[1:-1][0]:.2f}'))

    def plot_columns(df, x_col, y_col, title=None, xlabel=None, ylabel=None):
        """
        Plots two columns from a DataFrame against each other.

        Parameters:
            df (pd.DataFrame): The input DataFrame.
            x_col (str): The name of the column to use as the x-axis.
            y_col (str): The name of the column to use as the y-axis.
            title (str, optional): Title of the plot.
            xlabel (str, optional): Label for the x-axis. Defaults to x_col.
            ylabel (str, optional): Label for the y-axis. Defaults to y_col.
        """
        plt.figure(figsize=(10, 6))
        plt.plot(df[x_col], df[y_col], marker='o', linestyle='-', markersize=2)
        plt.title(title if title else f"{y_col} vs {x_col}")
        plt.xlabel(xlabel if xlabel else x_col)
        plt.ylabel(ylabel if ylabel else y_col)
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        

class universe_coord_extraction():
    # This class can extract coordinates and atom types of from mdanalysis universes and can write them to xyz files
    # Optional ability to make orca input files too

    def __init__(self):
        pass
 

    @staticmethod
    def mdanalysis_info_2_xyz(atom_types, coordinates, filepath):
        """
        Converts atom types and coordinates into an XYZ file format and saves it.
    
        Parameters:
            atom_types (list of str): List of atom type strings (e.g., 'h', 'c', 'o').
            coordinates (list of lists): List of coordinates, where each coordinate is [x, y, z].
            filepath (str): Path to save the XYZ file.
        """
        lines = []
    
        # Number of atoms
        no_atoms = len(atom_types)
        lines.append(f"{no_atoms}\n")  # First line: Number of atoms
        lines.append("\n")  # Second line: Blank comment line
    
        # Map lowercase atom types to proper case
        atom_map = {"h": "H", "c": "C", "o": "O"}
    
        # Generate the lines for each atom
        for i in range(len(atom_types)):
            atom_type = atom_types[i][0].lower()  # Convert to lowercase for mapping
            atom_name = atom_map.get(atom_type, atom_type.upper())  # Fallback to uppercase if not in map
            x, y, z = coordinates[i]  # Unpack coordinates
            line = f"{atom_name:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n"  # Format the line
            lines.append(line)
    
        # Write to the file
        with open(filepath, 'w') as file:
            file.writelines(lines)

        #print(f"XYZ file written to {filepath}")
        return(filepath)

    @staticmethod
    def save_polymers_xyz(universe_object, gen_orca_input=False):
        for i in range(len(universe_object.masterclass.poly_sel_dict)):
            polymer_name = "Polymer_" + str(i+1)
            poly = universe_object.select_polymer(polymer_name)
            coords, atom_types = poly.positions, poly.atoms.types

            xyz_filepath = os.path.join(universe_object.masterclass.poly_dft_input_dir, universe_object.masterclass.base_molecule_name + "_" + str(i+1) +  ".xyz")
            universe_coord_extraction.mdanalysis_info_2_xyz(atom_types, coords, xyz_filepath)
            if gen_orca_input == True:
                from modules.sw_file_formatter import DFT_input_generator
                input_filepath = os.path.join(universe_object.masterclass.poly_dft_input_dir, universe_object.masterclass.base_molecule_name + "_" + str(i+1) + ".inp")
                DFT_input_generator.generate_input(xyz_filepath, input_filepath, universe_object.masterclass.base_molecule_name + "_" + str(i+1))
               # print(f"Orca input file written to {input_filepath}")
        return(None)































