from modules.sw_orca import *
from modules.sw_basic_functions import *
import os as os
import subprocess
import csv
import math
import pandas as pd

class complex_fluid_model:
    """
    A standardized representation of a complex fluid model composed of multiple molecules.

    This class aggregates molecular data into a unified model format, computing weighted 
    averages of several properties and estimating the minimum requirements for simulation.

    This class is the output of a staticmethod in the 'complex_fluid_models' class.

    Attributes:
        molecules (list): A list of individual molecule objects.
        group_molecules (list): Grouped molecules, if applicable.
        molecule_ratios (list): Relative ratios of each molecule in the model.
        model_type (str): A string identifier for the type of model (e.g., "FT_model").
        model_name (str): A unique name for the model. Final name includes model_type as a suffix.
        wa_mw (float): Weighted average molecular weight.
        wa_chemical_hardness (float): Weighted average chemical hardness.
        wa_polarizability (float): Weighted average polarizability.
        wa_dipole_moment (float): Weighted average dipole moment.
        wa_total_energy (float): Weighted average total energy.
        wa_oxygen_content (float): Percentage of oxygen atoms across molecules.
        wa_nitrogen_content (float): Percentage of nitrogen atoms across molecules.
        wa_sulfur_content (float): Percentage of sulfur atoms across molecules.
        min_mols_for_sim (int): Minimum number of molecules required for a simulation.
        min_atoms_for_sim (int): Minimum number of atoms required for a simulation.
        min_vol_for_sim (float): Minimum volume required for a simulation.

    Args:
        molecules (list): List of molecule objects to include in the model.
        group_molecules (list): Grouped molecule objects (optional).
        molecule_ratios (list): List of floats indicating the contribution of each molecule.
        model_type (str): Type of the model used (e.g., "FT_model", "Cluster_model").
        model_name (str): Descriptive name for the model, will be suffixed with model_type.
    """
    def __init__(self, molecules, group_molecules, molecule_ratios, model_type, model_name):
        self.molecules = molecules
        self.group_molecules = group_molecules
        self.molecule_ratios = molecule_ratios
        self.model_type = model_type
        self.model_name = model_name
        if self.model_name != None:   
            self.model_name = model_name + "_" + model_type
        self.wa_mw = complex_fluid_models.get_weighted_average(self, "mw")
        self.wa_chemical_hardness = complex_fluid_models.get_weighted_average(self, "chemical_hardness")
        self.wa_polarizability = complex_fluid_models.get_weighted_average(self, "polarizability")
        self.wa_dipole_moment = complex_fluid_models.get_weighted_average(self, "dipole_moment")
        self.wa_total_energy = complex_fluid_models.get_weighted_average(self, "total_energy")
        self.wa_oxygen_content = complex_fluid_models.calculate_heteroatom_percentages(self.molecules)['O']
        self.wa_nitrogen_content = complex_fluid_models.calculate_heteroatom_percentages(self.molecules)['N']
        self.wa_sulfur_content = complex_fluid_models.calculate_heteroatom_percentages(self.molecules)['S']
        self.min_mols_for_sim = complex_fluid_models.min_mols_4_simulation(self)
        self.min_atoms_for_sim = complex_fluid_models.min_atoms_4_simulation(self)
        self.min_vol_for_sim = complex_fluid_models.min_vol_4_simulation(self)
        
class complex_fluid_models:
    """
    A class for generating and manipulating complex fluid models for chemical systems, as described in the referenced paper.

    This class contains a suite of functions designed to process and model complex fluid systems based on input 
    data from Orca molecules. The primary purpose of the class is to facilitate the generation of various 
    types of fluid models, each of which is standardized and ready for further analysis. The models generated 
    by this class are instances of the 'complex_fluid_model' class and represent different types of fluid systems 
    based on molecular properties, ratios, and interactions.

    The models include (but are not limited to):
    - Abundancy Grouped Model (AG_model)
    - Proportional Threshold Model (PT_model)
    - Scored Grouped Model (SG_model)
    - Others based on varying levels of molecular interactions and properties

    Key operations and features include:
    - Handling and processing Orca molecule data (e.g., peak areas, molecular weights, chemical hardness, etc.)
    - Grouping molecules based on various criteria (e.g., molecular similarity, peak area)
    - Calculating and normalizing molecular properties such as dipole moments, polarizability, and total energy
    - Generating various models based on these properties (e.g., AG_model, SG_model)
    - Providing a standardized output for further analysis or visualization

    The class is built to interface with the Orca molecular data and simplifies the process of modeling complex fluid systems.

    Methods include:
    - Generating different fluid models (e.g., AG_model, SG_model)
    - Performing calculations on molecular properties (e.g., weighted averages)
    - Grouping molecules and assigning ratios based on properties
    - Evaluating and comparing models using benchmark values
    - Generating dataframes with model results and rankings

    Attributes:
    -----------
    None (This class serves as a container for various static methods and does not store state itself)
    """
    def __init__(self):
        pass

    @staticmethod
    def get_weighted_average(model, class_attribute):
        """
        Calculate the weighted average of a specified molecular attribute across all molecules in the model.

        This function handles both:
        - Numerical attributes directly available on the molecule objects (e.g., "mw", "dipole_moment").
        - Heteroatom content (e.g., "Heteroatom_content_O") computed from SMILES strings.

        Args:
            model (complex_fluid_model): The complex fluid model instance containing molecules and their ratios.
            class_attribute (str): The name of the attribute to average.
                Use exact attribute name for molecule properties (e.g., "mw", "dipole_moment"),
                or one of the following for heteroatom content:
                - "Heteroatom_content_O"
                - "Heteroatom_content_N"
                - "Heteroatom_content_S"

        Returns:
            float: The weighted average value of the given attribute across all molecules.

        Raises:
            AttributeError: If the attribute doesn't exist on the molecule and isn't a supported heteroatom type.
        """
        if class_attribute.startswith("Heteroatom_content_"):
            atom_type = class_attribute.split("_")[-1].upper()
            atomic_weights = {"O": 16, "N": 14, "S": 32}

            if atom_type not in atomic_weights:
                raise ValueError(f"Unsupported heteroatom type: {atom_type}")

            def get_heteroatom_weight(smiles):
                return smiles.lower().count(atom_type.lower()) * atomic_weights[atom_type]

            weighted_sum = sum(
                ratio * get_heteroatom_weight(molecule.smiles)
                for ratio, molecule in zip(model.molecule_ratios, model.molecules)
            )
            return weighted_sum

        else:
            try:
                weighted_average = sum(
                    ratio * float(getattr(molecule, class_attribute))
                    for ratio, molecule in zip(model.molecule_ratios, model.molecules)
                )
                return weighted_average
            except AttributeError as e:
                raise AttributeError(
                    f"Attribute '{class_attribute}' not found in molecule objects."
                ) from e
            
    @staticmethod
    def calculate_heteroatom_percentages_sing_mol(molecule):
        """
        Calculate the percentage by weight of oxygen (O), nitrogen (N), and sulfur (S)
        atoms across a list of molecules, based on their SMILES strings.

        Args:
            molecule_list (list): List of molecule objects. Each must have a `.smiles` attribute.

        Returns:
            dict: A dictionary with keys 'O', 'N', and 'S' representing the total
                  mass percentage of each heteroatom across all molecules.

                  Example:
                  {
                      'O': 8.2,
                      'N': 5.1,
                      'S': 3.7
                  }

        Notes:
            - Uses atomic weights: O = 16, N = 14, S = 32.
            - SMILES strings are parsed naively (by counting letters).
            - This is an approximate method and does not handle full parsing of molecules.

        Raises:
            AttributeError: If a molecule in the list does not have a 'smiles' attribute.
        """
        smiles = molecule.smiles
        weight = molecule.mw
        atomic_weights = {'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'P': 30.97, 'S': 32.07, 'F': 18.99, 'Cl': 35.45, 'Br': 79.90, 'I': 126.90}
        heteroatoms = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']  # List of heteroatoms
    
        atom_weights = {}  
        for atom in heteroatoms:
            content = 0.0    
            count = smiles.count(atom)
            content += (count * atomic_weights[atom])/weight
            atom_weights[atom] = content
        
        tot_h_wt = 0.0
        mol = Chem.MolFromSmiles(smiles)
        mol_with_hydrogens = Chem.AddHs(mol)
        count = 0
        for atom in mol_with_hydrogens.GetAtoms():
                if atom.GetSymbol() == 'H':
                    count += 1
        content = (count * 1.0)/weight
        tot_h_wt += content
        atom_weights['H'] = tot_h_wt
    
        for key in atom_weights:
            if isinstance(atom_weights[key], (int, float)):
                atom_weights[key] *= 100
    
        return atom_weights
        
    @staticmethod
    def calculate_heteroatom_percentages(molecules):
        """
        Calculate the mass percentage of selected heteroatoms (including hydrogen) 
        in a mixture of molecules based on their SMILES strings and peak areas.

        This function:
        - Computes weighted contributions of each heteroatom by atomic mass and SMILES frequency.
        - Uses RDKit to explicitly count hydrogen atoms.
        - Normalizes all heteroatom weights to return percentage contributions.

        Args:
            molecules (list): A list of molecule objects.
                              Each must have `.smiles` (str) and `.peak_area` (float) attributes.

        Returns:
            dict: A dictionary mapping each atom ('C', 'H', 'O', etc.) to its
                  percentage contribution by mass in the total mixture.

                  Example:
                  {
                      'C': 68.2,
                      'H': 10.3,
                      'O': 15.0,
                      'N': 6.5,
                      ...
                  }

        Requires:
            - RDKit must be installed and imported as `from rdkit import Chem`.

        Raises:
            AttributeError: If molecule objects do not have 'smiles' or 'peak_area' attributes.
        """
        # Extract SMILES and peak area for each molecule
        try:
            mixture = {mol.smiles: float(mol.peak_area) for mol in molecules}
        except AttributeError as e:
            raise AttributeError("Each molecule must have 'smiles' and 'peak_area' attributes.") from e
            
        total_peak_area = complex_fluid_models.get_group_area(molecules)
    
        atomic_weights = {'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'P': 30.97, 'S': 32.07, 'F': 18.99, 'Cl': 35.45, 'Br': 79.90, 'I': 126.90}
        heteroatoms = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']  # List of heteroatoms
      
        atom_weights = {}  
        for atom in heteroatoms:
            content = 0.0
            for compound, weight in mixture.items():
                count = compound.count(atom)
                content += (weight/total_peak_area) * count * atomic_weights[atom]
            atom_weights[atom] = content
         # Now need to get H content
        tot_h_wt = 0.0
        for compound, weight in mixture.items():
            hydrogen_weight = 0.0
            mol = Chem.MolFromSmiles(compound)
            mol_with_hydrogens = Chem.AddHs(mol)
            count = 0
            for atom in mol_with_hydrogens.GetAtoms():
                    if atom.GetSymbol() == 'H':
                        count += 1
            content = (weight/total_peak_area) * count * 1.0
            tot_h_wt += content
        atom_weights['H'] = tot_h_wt
        heteroatoms.append('H') 
        total_weight = sum(atom_weights.values())
        for atom in heteroatoms:
            if total_weight == 0:
                atom_weights[atom] = 0
            else:
                atom_weights[atom] = (atom_weights[atom]/total_weight)*100
   
        return atom_weights # this returns the heteroatom percentage 
    
    @staticmethod
    def group_molecules(molecule_class_list):
        """
        Groups molecules into categories based on the presence of nitrogen and/or oxygen 
        heteroatoms and the types of rings (5-membered, 6-membered, mixed, or no ring) 
        present in their structures.

        Molecules are categorized into 4 main groups:
        1. Contains both nitrogen and oxygen
        2. Contains only nitrogen
        3. Contains only oxygen
        4. Contains neither nitrogen nor oxygen

        Within each group, molecules are further divided into:
        - [0] 5-membered ring only
        - [1] 6-membered ring only
        - [2] mixed or other rings
        - [3] no rings

        Args:
            molecule_class_list (list): List of molecule objects with a `.smiles` attribute.

        Returns:
            list: A nested list structure:
                  [mols_with_n_and_o, mols_with_n, mols_with_o, no_hetero]
                  where each sublist contains 4 subcategories as described above.

        Requires:
            - `has_heteroatoms(mol)` to return a list of heteroatoms present.
            - `has_rings(mol)` to return a list of detected ring types.
            - RDKit to be installed for SMILES parsing via `Chem.MolFromSmiles`.
        """
        smiles_list = [molecule.smiles for molecule in molecule_class_list]
        
        mols_with_n_and_o = [[], [], [], []]       
        mols_with_n = [[], [], [], []] # 5 RING, 6 RING, OTHER RING(MAYBE COMBO CORES), NO RING
        mols_with_o = [[], [], [], []]
        no_hetero = [[], [], [], []]      
        N = "nitrogen"
        O = "oxygen"
        S = "sulfur"

        five = "5-membered ring"
        six = "6-membered ring"
     
        for molecule in molecule_class_list:
            mol = Chem.MolFromSmiles(molecule.smiles)
        
            a = has_heteroatoms(mol)
            if N in a:
                if O in a:
                    b = has_rings(mol)
                    if five in b and six not in b:
                        mols_with_n_and_o[0].append(molecule) # Pulls out 5 membered rings
                    elif six in b and five not in b:
                        mols_with_n_and_o[1].append(molecule) # Pulls out 6 membered rings
                    elif len(b) != 0:
                        mols_with_n_and_o[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
                    else:
                        mols_with_n_and_o[3].append(molecule)
                elif O not in a:
                    b = has_rings(mol)
                    if five in b and six not in b:
                        mols_with_n[0].append(molecule) # Pulls out 5 membered rings
                    elif six in b and five not in b:
                        mols_with_n[1].append(molecule) # Pulls out 6 membered rings
                    elif len(b) != 0:
                        mols_with_n[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
                    else:
                        mols_with_n[3].append(molecule)           
                else:
                    pass
            elif O in a:
                b = has_rings(mol)
                if five in b and six not in b:
                    mols_with_o[0].append(molecule) # Pulls out 5 membered rings
                elif six in b and five not in b:
                    mols_with_o[1].append(molecule) # Pulls out 6 membered rings
                elif len(b) != 0:
                    mols_with_o[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
                else:
                    mols_with_o[3].append(molecule)
            else:
                b = has_rings(mol)
                if five in b and six not in b:
                    no_hetero[0].append(molecule) # Pulls out 5 membered rings
                elif six in b and five not in b:
                    no_hetero[1].append(molecule) # Pulls out 6 membered rings
                elif len(b) != 0:
                    no_hetero[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
                else:
                    no_hetero[3].append(molecule)
        grouped_mols = []
        grouped_mols.append(mols_with_n_and_o)
        grouped_mols.append(mols_with_n)
        grouped_mols.append(mols_with_o)
        grouped_mols.append(no_hetero)
        return(grouped_mols)
    
    @staticmethod
    def is_orca_molecule(molecule):
        """
        Determines whether the provided object is an instance of the `orca_molecule` class.

        This utility method helps validate that a given molecule belongs to the expected 
        class used for ORCA-related molecular data.

        Args:
            molecule (object): The object to check.

        Returns:
            bool: True if the object is an instance of `orca_molecule`, False otherwise.
        """
        return isinstance(molecule, orca_molecule)

    @staticmethod
    def all_model(model_name, orca_molecules):
        """
        Creates a `complex_fluid_model` using all provided `orca_molecule` instances with 
        their peak areas used to determine composition ratios.

        This model includes all molecules regardless of thresholds or grouping, and assumes
        that the `peak_area` values are percentages that sum to ~100.

        Args:
            model_name (str): A name for the model being generated.
            orca_molecules (list): A list of `orca_molecule` instances to include in the model.

        Returns:
            complex_fluid_model or None: Returns a `complex_fluid_model` instance if all inputs 
            are valid `orca_molecule` objects. If validation fails, prints an error and returns `None`.
        """

        # Validate input: all items must be instances of orca_molecule
        if not all(isinstance(molecule, orca_molecule) for molecule in orca_molecules):
            print("Please use orca_molecule class instances for each molecule in the list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return None

        # Calculate ratios from peak areas assuming they are given as percentages
        ratios = [float(mol.peak_area) / 100 for mol in orca_molecules]

        orca_mol_bool = [isinstance(molecule, orca_molecule) for molecule in orca_molecules]
      
        # Calculate ratios from peak areas assuming they are given as percentages
        ratios = [float(mol.peak_area) / 100 for mol in orca_molecules]
        all_model = complex_fluid_model(orca_molecules, None, ratios, "ALL_model", model_name)
        
        return(all_model)

    @staticmethod
    def fixed_threshold_model(model_name=None, orca_molecules=None, selection_threshold=None):
        """
        Creates a `complex_fluid_model` by filtering `orca_molecule` instances based on a peak area threshold.

        This method filters out any molecules whose `peak_area` is below the specified threshold and 
        constructs a model only from those that meet or exceed it.

        Args:
            model_name (str): A name for the resulting model.
            orca_molecules (list): A list of `orca_molecule` instances.
            selection_threshold (float): The minimum peak area required for a molecule to be included.

        Returns:
            complex_fluid_model or None: A `complex_fluid_model` instance created from filtered molecules. 
            Returns `None` if input validation fails.
        """
        # Input validation
        if None in (model_name, orca_molecules, selection_threshold):
            print("Error: Please ensure all arguments are provided.")
            print("EXAMPLE: ft_model = complex_fluid_models.fixed_threshold_model("
              "model_name='model_name', orca_molecules=list_of_molecule_objects, selection_threshold=10)")
            return None

        # Check all objects are instances of orca_molecule
        if not all(isinstance(molecule, orca_molecule) for molecule in orca_molecules):
            print("Please use orca_molecule class instances for each molecule in the list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return None

        fixed_threshold_model = [] 

        for molecule in orca_molecules:
            if float(molecule.peak_area) >= selection_threshold:
                fixed_threshold_model.append(molecule)

        if not fixed_threshold_model:
            print("No molecules met the peak area threshold. Model not generated.")
            return(None)
        
        tot_peak_area = sum([float(mol.peak_area) for mol in fixed_threshold_model])
        ratios = [(float(mol.peak_area)/tot_peak_area) for mol in fixed_threshold_model]
        
        fixed_threshold_model = complex_fluid_model(fixed_threshold_model, None, ratios, "FT_model", model_name)
        
        return(fixed_threshold_model)

    @staticmethod
    def proportional_threshold_model(model_name=None, orca_molecules=None):
        """
        Constructs a `complex_fluid_model` by selecting a proportionally normalized set of molecules.

        This method ensures that the number of selected molecules does not exceed the original number
        of molecules provided, using a normalization strategy based on the ratio of each molecule's 
        peak area to the minimum peak area in the list.

        Args:
            model_name (str): Name to assign to the resulting model.
            orca_molecules (list): A list of `orca_molecule` instances.

        Returns:
            complex_fluid_model or None: A `complex_fluid_model` instance composed of proportionally
            selected molecules, or `None` if validation fails.
        """

        # Input validation
        if None in (model_name, orca_molecules):
            print("Error: please ensure all arguments are passed")
            print("EXAMPLE: pt_model = complex_fluid_models.proportional_threshold_model("
                  "model_name='model_name', orca_molecules=list_of_molecule_objects)")
            return None

        # Type validation
        if not all(isinstance(molecule, orca_molecule) for molecule in orca_molecules):
            print("Please use orca_molecule class instances for each molecule in list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return None

        def calculate_normalised_total(all_info):
            normalised_values = []
            for i in range(len(all_info)):
                normalised_values.append(all_info[i][3])
            return(sum(normalised_values))
        
        all_info = []
        peak_areas = [] # Used for initial normalising calculation

        for molecule in orca_molecules:
            peak_areas.append(float(molecule.peak_area))

        max_mols = len(orca_molecules) # The number of molecules to select for the model
        
        for molecule in orca_molecules:
            info = [molecule, molecule.smiles, molecule.peak_area, (float(molecule.peak_area)/min(peak_areas))]
            all_info.append(info)

        for i in range(len(orca_molecules)):
            normalised_total = calculate_normalised_total(all_info)
            if normalised_total > max_mols:
                min_peak_area = peak_areas.index(min(peak_areas))
                peak_areas.pop(min_peak_area)
                all_info.pop(min_peak_area)

            for j in range(len(all_info)):
                all_info[j][3] = float(all_info[j][2])/min(peak_areas)
  
        selected_mols = []
        
        for i in range(len(all_info)):
            selected_mols.append(all_info[i][0])

        tot_peak_area = sum([float(mol.peak_area) for mol in selected_mols])
        ratios = [(float(mol.peak_area)/tot_peak_area) for mol in selected_mols]

        proportional_threshold_model = complex_fluid_model(selected_mols, None, ratios, "PT_model", model_name)

        return(proportional_threshold_model)

    @staticmethod
    def get_group_area(orca_molecules):
        """
        Calculates the total peak area for a list of `orca_molecule` instances.

        This method sums the `peak_area` attribute of each molecule in the provided list
        and returns the total peak area, which can be used for further calculations in
        models like AG and SG.

        Args:
            orca_molecules (list): A list of `orca_molecule` instances. Each molecule
            is expected to have a `peak_area` attribute.

        Returns:
            float: The total peak area computed as the sum of the `peak_area` attributes
            of all molecules in the input list.
        """
        total_peak_area = sum(float(molecule.peak_area) for molecule in orca_molecules)
        return total_peak_area
        
    @staticmethod
    def abundancy_grouped_model(model_name=None, orca_molecules=None):
        """
        Generates an abundance-based grouped model from a list of `orca_molecule` instances.

        This method groups molecules based on their structural features, selects the
        molecule with the largest peak area from each group, and generates a model using 
        these representative molecules. The relative abundance of each group is computed 
        based on the total peak area of molecules in that group.

        Args:
            model_name (str, optional): The name of the model to be created. If not provided, 
            the model will not be named.
            orca_molecules (list): A list of `orca_molecule` instances. Each molecule is
            expected to have a `peak_area` attribute.

        Returns:
            complex_fluid_model: A `complex_fluid_model` instance created using the selected
            representative molecules, their ratios, and other related data.

        Raises:
            ValueError: If any of the arguments are `None` or if invalid `orca_molecule` 
            instances are provided.
        """
        # Check if model_name and orca_molecules are provided
        if None in (model_name, orca_molecules):
            print("Error: please ensure all arguments are passed")
            print("EXAMPLE: ag_model = complex_fluid_models.abundancy_grouped_model(model_name='model_name', orca_molecules=list_of_molecule_objects)")
            return None  

        # Validate that each molecule in the list is an instance of orca_molecule
        orca_mol_bool = [isinstance(molecule, orca_molecule) for molecule in orca_molecules]
        if "False" in orca_mol_bool:
            print("Please use orca molecule class instances for each molecule in list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return None

        grouped_molecules = complex_fluid_models.group_molecules(orca_molecules)

        abundancy_grouped_model = []
        all_mols_in_group = []

        groups_in_model = 0 
        for j in range(len(grouped_molecules)):
            for k in range(len(grouped_molecules[j])):
                group = grouped_molecules[j][k]
                if len(group) == 0:
                    pass
                else:
                    all_mols_in_group.append(group)
                    groups_in_model += 1
                    molecule_to_model = group[0]
                    for i in range(len(group)):
                        if i == 0:
                            pass
                        else:
                            if float(group[i].peak_area) > float(molecule_to_model.peak_area):
                                molecule_to_model = group[i]
                            else:
                                pass
                    abundancy_grouped_model.append(molecule_to_model)

        group_peaks = []
        mw = []

        for j in range(len(grouped_molecules)):
            for k in range(len(grouped_molecules[j])):
                group = grouped_molecules[j][k]
                if len(group) == 0:
                    pass
                else:
                    group_area = complex_fluid_models.get_group_area(group)
                    group_peaks.append(group_area)

        summ = sum(group_peaks)
        for i in range(len(group_peaks)):
            group_peaks[i] = group_peaks[i]/summ

        abundancy_grouped_model = complex_fluid_model(abundancy_grouped_model, all_mols_in_group, group_peaks, "AG_model", model_name)

        return(abundancy_grouped_model)
        
    @staticmethod
    def scored_grouped_model(model_name=None, orca_molecules=None):
        """
        Generates a scored grouped model based on the selection of the best molecules
        from grouped molecule categories, with scores calculated based on a variety of
        molecular properties.

        This method uses the concept of an 'abundance-grouped model' (AG model) as a base,
        and then applies scoring to each molecule in each group. The molecule with the best
        score is selected from each group, where the score is determined by the weighted
        differences between normalized molecular properties like molecular weight, total
        energy, polarizability, dipole moment, chemical hardness, and oxygen content.

        Args:
            model_name (str, optional): The name of the model to be created. If not provided, 
            the model will not be named.
            orca_molecules (list): A list of `orca_molecule` instances. Each molecule is
            expected to have the required attributes such as `peak_area`, `mw`, `total_energy`,
            `polarizability`, `dipole_moment`, `chemical_hardness`, and methods to calculate 
            heteroatom percentages.

        Returns:
            complex_fluid_model: A `complex_fluid_model` instance created using the selected
            representative molecules, their ratios, and other related data.

        Raises:
            ValueError: If any of the arguments are `None` or if invalid `orca_molecule` 
            instances are provided.
        """
        # Check if model_name and orca_molecules are provided
        if None in (model_name, orca_molecules):
            print("Error: please ensure all arguments are passed")
            print("EXAMPLE: sg_model = complex_fluid_models.scored_grouped_model(model_name='model_name', orca_molecules=list_of_molecule_objects)")
            return None  

        
        all_model = complex_fluid_models.all_model(model_name, orca_molecules)
        # this is a variant of the ag_model
        ag_model = complex_fluid_models.abundancy_grouped_model(model_name, orca_molecules)

        ag_model = ag_model.molecules, ag_model.group_molecules, ag_model.molecule_ratios 

        selected_molecules = []
        all_groups = []
        group_peaks = []

        for i in range(len(ag_model[1])):
            scored_molecules = [[], []]
            model = complex_fluid_models.all_model(None, ag_model[1][i])
            print(model)
            all_groups.append(model)
            group_peak = complex_fluid_models.get_group_area(model.molecules)
            group_peaks.append(group_peak)

            if len(model.molecules) == 1:
                scored_molecules[0].append(model.molecules[0])
                scored_molecules[1].append(1.0)
                selected_molecules.append(model.molecules[0])
                continue
            print(model)
            # These don't need to be calculated here - can be called from the instance
            oxygen_content = complex_fluid_models.calculate_heteroatom_percentages(model.molecules)['O']
            wa_mw = complex_fluid_models.get_weighted_average(model, "mw")
            wa_tot_en = complex_fluid_models.get_weighted_average(model, "total_energy")
            wa_polar = complex_fluid_models.get_weighted_average(model, "polarizability")
            wa_dipole = complex_fluid_models.get_weighted_average(model, "dipole_moment")
            wa_chem_hard = complex_fluid_models.get_weighted_average(model, "chemical_hardness")

            mws = [float(mol.mw) for mol in model.molecules]
            min_mw = min(mws)
            max_mw = max(mws)

            normalized_mws = [(mw - min_mw) / (max_mw - min_mw) for mw in mws]
            avg_norm_mw = (wa_mw - min_mw) / (max_mw - min_mw)
        
            tot_ens = [float(mol.total_energy) for mol in model.molecules]
            min_tot_en = min(tot_ens)
            max_tot_en = max(tot_ens)
            normalized_tot_ens = [(tot_en - min_tot_en) / (max_tot_en - min_tot_en) for tot_en in tot_ens]
            avg_norm_tot_en = (wa_tot_en - min_tot_en) / (max_tot_en - min_tot_en)
        
            polars = [float(mol.polarizability) for mol in model.molecules]
            min_polar = min(polars)
            max_polar = max(polars)
            normalized_polars = [(polar - min_polar) / (max_polar - min_polar) for polar in polars]
            avg_norm_polar = (wa_polar - min_polar) / (max_polar - min_polar)        
                                            
            dipoles = [float(mol.dipole_moment) for mol in model.molecules]
            min_dipole = min(dipoles)
            max_dipole = max(dipoles)
            normalized_dipoles = [(dipole - min_dipole) / (max_dipole - min_dipole) for dipole in dipoles]
            avg_norm_dipole = (wa_dipole - min_dipole) / (max_dipole - min_dipole)
                                            
            chem_hards = [float(mol.chemical_hardness) for mol in model.molecules]
            min_chem_hard = min(chem_hards)
            max_chem_hard = max(chem_hards)
            normalized_chem_hards = [(chem_hard - min_chem_hard) / (max_chem_hard - min_chem_hard) for chem_hard in chem_hards]
            avg_norm_chem_hard = (wa_chem_hard - min_chem_hard) / (max_chem_hard - min_chem_hard)
        
            oxys = [complex_fluid_models.calculate_heteroatom_percentages_sing_mol(mol)['O'] for mol in model.molecules]
            min_oxy = min(oxys)
            max_oxy = max(oxys)
            if min_oxy == 0:
                if max_oxy ==0:
                    key="no_oxy"           
            else:
                key="oxy_present"
                normalized_oxys = [(oxy - min_oxy ) / (max_oxy  - min_oxy ) for oxy  in oxys]
                avg_norm_oxy = (oxygen_content - min_oxy) / (max_oxy - min_oxy)
        
            for j in range(len(model.molecules)):
                # This key on;y exists to to create an alternative calcualtion for those molecules that do not have oxygens
                if key == "oxy_present":
                    score_value = 1/(1+(math.sqrt((normalized_mws[j]-float(avg_norm_mw))**2 + (normalized_tot_ens[j]-float(avg_norm_tot_en))**2 +(normalized_polars[j]-float(avg_norm_polar))**2 + (normalized_dipoles[j]-float(avg_norm_dipole))**2 +(normalized_chem_hards[j]-float(avg_norm_chem_hard))**2 +(normalized_oxys[j]-float(avg_norm_oxy))**2)))
                else:
                    score_value = 1/(1+(math.sqrt((normalized_mws[j]-float(avg_norm_mw))**2 + (normalized_tot_ens[j]-float(avg_norm_tot_en))**2 +(normalized_polars[j]-float(avg_norm_polar))**2 + (normalized_dipoles[j]-float(avg_norm_dipole))**2 +(normalized_chem_hards[j]-float(avg_norm_chem_hard))**2)))
            
                scored_molecules[0].append(model.molecules[j])
                scored_molecules[1].append(score_value)
        
            best_score = max(scored_molecules[1])
            best_score_index = scored_molecules[1].index(best_score)

            selected_molecules.append(scored_molecules[0][best_score_index])  
        group_peaks = [peak/100 for peak in group_peaks]
        scored_grouped_model = complex_fluid_model(selected_molecules, all_groups, group_peaks, "SG_model", model_name)
        
        return(scored_grouped_model)

    @staticmethod
    def generate_model_df(models): # molecules is all set, models are subsets of molecules
        """
        Generates a DataFrame containing various properties for a list of molecular models.

        This method extracts key properties such as molecular weight (Mw), chemical hardness,
        polarizability, dipole moment, total energy, and oxygen content for each model in the 
        provided list. These properties are calculated as weighted averages over the molecules 
        within each model. The resulting DataFrame organizes these properties along with the 
        model type for each model.

        Args:
            models (list): A list of `complex_fluid_model` instances, where each model contains
                           a set of molecules and is associated with a particular model type. 
                           Each model should have molecules and associated properties 
                           (like `model_type` and `molecules`).

        Returns:
            pd.DataFrame: A pandas DataFrame containing the following columns:
                - "Model Type": The type of each model.
                - "Mw": The weighted average molecular weight for each model.
                - "Chem_hard": The weighted average chemical hardness for each model.
                - "polarizability": The weighted average polarizability for each model.
                - "Dipole": The weighted average dipole moment for each model.
                - "total_energy": The weighted average total energy for each model.
                - "oxygen_content": The weighted average oxygen content for each model.
    
        Raises:
            ValueError: If any of the models in the list do not have the required attributes 
                            or are not valid `complex_fluid_model` instances.
        """
        mws = []
        chem_hards = []
        polars = []
        dipoles = []
        tot_ens = []
        oxygen_content = []
        model_types = []
    
        for i in range(len(models)):
            mws.append(complex_fluid_models.get_weighted_average(models[i], "mw"))
            chem_hards.append(complex_fluid_models.get_weighted_average(models[i], "chemical_hardness"))
            polars.append(complex_fluid_models.get_weighted_average(models[i], "polarizability"))
            dipoles.append(complex_fluid_models.get_weighted_average(models[i], "dipole_moment"))
            tot_ens.append(complex_fluid_models.get_weighted_average(models[i], "total_energy"))
            model_types.append(models[i].model_type)

            oxy = complex_fluid_models.calculate_heteroatom_percentages(models[i].molecules)["O"]
            oxygen_content.append(oxy)
    
        # Create a DataFrame with model properties
        df = pd.DataFrame({
            "Model Type": model_types,
            "Mw": mws,
            "Chem_hard": chem_hards,
            "polarizability": polars,
            "Dipole": dipoles,
            "total_energy": tot_ens,
            "oxygen_content": oxygen_content
        })
    
        return df

    @staticmethod
    def rank_models_dep(df):
        # Ensure the column names are correct and in proper case
        df.columns = df.columns.str.strip()
    
        properties = ["Mw", "Chem_hard", "polarizability", "Dipole", "total_energy", "oxygen_content"]

        # Create a new dataframe to store the calculated differences and ranks
        differences_df = pd.DataFrame()
        ranks_df = pd.DataFrame()
    
        # Copy the "Model Type" column from the input dataframe to the output dataframes
        differences_df["Model Type"] = df["Model Type"]
        ranks_df["Model Type"] = df["Model Type"]

        # Check if 'ALL_model' is in the "Model Type" column
        if "ALL_model" not in df["Model Type"].values:
            raise ValueError("ALL_model not found in 'Model Type' column")

        # Calculate absolute differences for each property against the "all_model" row
        for prop in properties:
            all_model_value = df.loc[df["Model Type"] == "ALL_model", prop].iloc[0]
            differences_df[prop + "_Difference"] = df[prop].sub(all_model_value).abs()

            # Calculate ranks for each descriptor against the "all_model" benchmark
            ranks_df[prop + "_Rank"] = differences_df[prop + "_Difference"].rank()

        # Calculate cumulative difference for each model based on the property differences
        differences_df["Cumulative_Difference"] = differences_df[[prop + "_Difference" for prop in properties]].sum(axis=1)

        # Calculate final score for each model
        min_diff = differences_df["Cumulative_Difference"].min()
        differences_df["Final_Score"] = differences_df["Cumulative_Difference"].apply(lambda diff: (min_diff - diff) + 1)

        # Calculate the sum of ranks for each model
        ranks_df["Final_Rank"] = ranks_df[[prop + "_Rank" for prop in properties]].sum(axis=1)

        return ranks_df

    @staticmethod
    def rank_models(model_df, benchmark_model_idx=0, features=None):
        """
        Calculates the Euclidean distance of each model from the benchmark model
        based on multiple features and ranks the models accordingly.
    
        Args:
            model_df (pd.DataFrame): The DataFrame containing model properties.
            benchmark_model_idx (int): The index of the benchmark model in the DataFrame.
            features (list): List of feature columns to compare (e.g., ['Mw', 'Chem_hard', 'polarizability']).

        Returns:
            pd.DataFrame: DataFrame with the models ranked based on their distance to the benchmark model.
        """
        if features is None:
            features = ["Mw", "Chem_hard", "polarizability", "Dipole", "total_energy", "oxygen_content"]
    
        # Get the feature values of the benchmark model
        benchmark_model = model_df.iloc[benchmark_model_idx][features]
    
        # Function to calculate Euclidean distance between model and benchmark model
        def euclidean_distance(model, benchmark):
            return math.sqrt(sum([(model[feature] - benchmark[feature])**2 for feature in features]))
    
        # Apply the distance calculation to each model in the DataFrame
        model_df['Distance_to_Benchmark'] = model_df.apply(lambda row: euclidean_distance(row, benchmark_model), axis=1)
    
        # Rank the models based on their distance (smallest distance is best match)
        model_df['Rank'] = model_df['Distance_to_Benchmark'].rank(ascending=True, method='min')
    
        # Sort by rank (best match first)
        ranked_df = model_df.sort_values(by='Rank', ascending=True)
    
        return ranked_df
        
    @staticmethod
    def min_mols_4_simulation(model):
        """
        Calculate the minimum number of molecules required for the simulation based on the 
        molecule ratios in the model.

        This method computes the minimum number of molecules (min_mols) required for a simulation 
        by evaluating the relative molecule ratios and adjusting them based on the total number of 
        molecules in the model. The number of molecules is scaled by the ratio of the molecule's 
        ratio to the minimum ratio in the model.

        Args:
            model (complex_fluid_model): The model object containing the list of molecules and their respective ratios.

        Returns:
            int: The minimum number of molecules (rounded down) required for the simulation.

        Example:
            min_mols = complex_fluid_models.min_mols_4_simulation(model)
        """
        num_compounds = len(model.molecules)
        min_peak = min(model.molecule_ratios)
        min_mols = 0
        for i in range(num_compounds):
            min_mols += (float(model.molecule_ratios[i]) / float(min_peak)) * (len(model.molecules) + 1)
    
        return int(min_mols)

    @staticmethod
    def min_atoms_4_simulation(model):
        """
        Calculate the minimum number of atoms required for the simulation based on the molecule ratios 
        and the molecular structure in the model.

        This method computes the minimum number of atoms required for the simulation by iterating 
        through each molecule in the model, retrieving its atomic structure, and adjusting the 
        number of atoms based on the ratio of the molecule's ratio to the minimum ratio in the model. 
        The number of atoms is then scaled by the total number of molecules in the model.

        Args:
            model (complex_fluid_model): The model object containing the list of molecules and their respective ratios.

        Returns:
            int: The minimum number of atoms (rounded down) required for the simulation.

        Example:
            min_atoms = complex_fluid_models.min_atoms_4_simulation(model)
        """
        min_peak = min(model.molecule_ratios)
        num_atoms = 0
        for i in range(len(model.molecules)):
            mol = Chem.MolFromSmiles(model.molecules[i].smiles)
            mol = Chem.AddHs(mol)
            num_atoms = num_atoms + (mol.GetNumAtoms()*(float(model.molecule_ratios[i]) / float(min_peak))*(len(model.molecules) + 1))
            
        return int(num_atoms)
        
    @staticmethod
    def min_vol_4_simulation(model):
        """
        Calculate the minimum total volume required for the simulation based on the molecule ratios and 
        the minimum number of molecules for the simulation.

        This method computes the total volume by multiplying each molecule's volume by its corresponding 
        ratio and the minimum number of molecules required for the simulation, and then summing these values 
        across all molecules in the model.

        Parameters:
        -----------
        model : object
            An instance of a model class containing the following attributes:
                - `molecules` (list): A list of molecule objects, each having a `volume` attribute.
                - `molecule_ratios` (list): A list of ratios corresponding to each molecule.
                - `min_mols_for_sim` (float): The minimum number of molecules required for the simulation.

        Returns:
        --------
        float
            The total minimum volume (in the same units as the molecule volumes) required for the simulation.
        """
        tot_volume = 0      
        for i in range(len(model.molecules)):
            volume_of_mol = model.molecules[i].volume*(model.molecule_ratios[i]*model.min_mols_for_sim)
            tot_volume = tot_volume + volume_of_mol
         
        return tot_volume

    @staticmethod
    def model_output_block(all_model, model, ranked_data):
        """
        Generates a detailed output block for a given model, comparing it with the benchmark model (ALL_model).
        The output includes key weighted values, heteroatom content, minimum molecule/atom requirements, and a
        representation of the molecules in the model.

        This function evaluates how well a specific model (e.g., FT_model, AG_model, PT_model) aligns with the 
        benchmark model (ALL_model). It calculates various properties (e.g., molecular weight, chemical hardness, 
        polarizability) for both the model and the benchmark, and includes additional details such as molecule 
        ratios, heteroatom content, and the minimum number of molecules and atoms required for unbiased simulation.

        Parameters:
        -----------
        all_model : object
            An instance of the "ALL_model", serving as the benchmark for comparison. The benchmark model's weighted 
            averages (e.g., molecular weight, chemical hardness) are used for scoring the provided model.
        
        model : object
            An instance of a model (e.g., FT_model, AG_model, PT_model) to be compared to the benchmark (ALL_model).
            This model contains various molecular properties and ratios for which comparison and evaluation are made.

        ranked_data : pandas.DataFrame
            A DataFrame containing the rankings of various models. This is used to get the score of the given model 
            against the benchmark for performance evaluation.

        Returns:
        --------
        list
            A list of strings representing a formatted output block detailing the model properties, including the 
            model's weighted values, heteroatom content, minimum number of molecules and atoms required for simulation,
            and the list of molecules in the model along with their SMILES notation.
        """
        if all_model.model_type != "ALL_model":
            print("Please provide the ALL_model to serve as a benchmark")
            return(None)
        supported_models = ["FT_model", "ALL_model", "AG_model", "PT_model", "SG_model"]
        if model.model_type not in supported_models:
            print("This function currently does not support: ", model.model_type)
            return(None)
        # Set up a dictionary for the benchmark so I can score each model versus the benchmark
        benchmark_values = {"Mw":0, "Chem_hard":0, "polarizability":0, "Dipole":0, "total_energy":0, "oxygen_content":0}
        benchmark_values["Mw"] = all_model.wa_mw
        benchmark_values["Chem_hard"] = all_model.wa_chemical_hardness
        benchmark_values["polarizability"] = all_model.wa_polarizability
        benchmark_values["Dipole"] = all_model.wa_dipole_moment
        benchmark_values["total_energy"] = all_model.wa_total_energy
        benchmark_values["oxygen_content"] = all_model.wa_oxygen_content
    
        # Creation of a dictionary for the model_values
        model_values = {"Mw":0, "Chem_hard":0, "polarizability":0, "Dipole":0, "total_energy":0, "oxygen_content":0} 
    
        num_mols = len(model.molecules)
          
        spacer = ""
        molecule_header = "The molecules in this model are:"
        smiles_header = "The smiles for these molecules are:"
        model_lines = []
        recorded_peak_area = 0
    
        # This loop gets total peak area of the entire dataset passed to it
        dataset_peak_area = sum(float(molecule.peak_area) for molecule in all_model.molecules)
        # This gets the peak area of the compounds in the model (or subset)
        model_area = sum(float(molecule.peak_area) for molecule in model.molecules)
    
        model_line = "Model: " + model.model_type
             
        # Representation of model in a string to write to file
        representation_line = "This model represents " + str((model_area/dataset_peak_area)*100) + "% of the charecterised biocrude"
        # Other properties of model to write to line
        weighted_line = "Weighted values for molecular weight, chemical hardness, polarizability, chemical hardness and dipole moment are listed below:"
        mw = "Molecular weight  = " + str(model.wa_mw)
        chem_hard = "Chemical hardness = " + str(model.wa_chemical_hardness)
        polar = "Polarizability    = " + str(model.wa_polarizability)
        dipole = "Dipole moment     = " + str(model.wa_dipole_moment)
        tot_energy = "Total energy      = " + str(model.wa_total_energy)
    
        model_values["oxygen_content"] = model.wa_oxygen_content            
        model_values["Mw"] = model.wa_mw
        model_values["Chem_hard"] = model.wa_chemical_hardness
        model_values["polarizability"] = model.wa_polarizability
        model_values["Dipole"] = model.wa_dipole_moment
        model_values["total_energy"] = model.wa_total_energy
    
        # Functions to calculate how good a model is versus the benchmark  
        model_score = ranked_data.loc[ranked_data["Model Type"] == model.model_type, "Final_Rank"].values[0]
        score_line = "This model has a score of: " + str(model_score/6) + " versus the benchmark"
        hetero_content_line = "The heteroatom content of this model is:"
        model_lines.extend([model_line, representation_line, score_line, spacer, weighted_line, mw, chem_hard, polar, dipole, tot_energy, spacer, hetero_content_line])
     
        # Now we write heteroatoms content to file
        heteroatom_content = complex_fluid_models.calculate_heteroatom_percentages(model.molecules)

        for atom, weight in heteroatom_content.items():
            atom_content = f"{atom}: {weight}"
            model_lines.append(atom_content)
        
        atoms_line = "This model will require the following minimum number of atoms for unbiased simulation:"
        mol_line = "This model will require the following minimum number of molecules for unbiased simulation:"
        model_lines.extend([spacer, atoms_line, str(model.min_atoms_for_sim), spacer, mol_line, str(model.min_mols_for_sim), spacer, molecule_header])
                           
        for i in range(num_mols):
            if model.model_type == "AG_model" or model.model_type == "SG_model":
                group_area = model.molecule_ratios[i]
                line = model.molecules[i].name + " and this group accounts for: " + str(group_area) + " % of the biocrude"
    
            else:
                model_lines.append(model.molecules[i].name)
           
        model_lines.extend([spacer, smiles_header])
        for i in range(num_mols):
            model_lines.append(model.molecules[i].smiles)
       
        return(model_lines) 

class complex_fluid_model_builder:
    # This  class contains all the functions to generate complex_fluid_models as described in the paper
    # The inputs are orca molecules
    # The outputs are instances of the 'complex_fluid_model' class - all standardised so it is easier to work with
    
    def __init__(self):
        pass

    @staticmethod
    def generate_packmol_bio_oil_cube(model_dirs, model, tolerance=None, filetype=None, volume_scalar=None, molecule_scalar=None): # csv_name (i.e. pine bark phenolic), molecules selected by model generation, model type
        total_volume = model.min_vol_for_sim
        
        volume_scalar = 1.5 if volume_scalar is None else volume_scalar

        molecule_scalar = 1 if molecule_scalar is None else molecule_scalar
        
        size = int((math.pow(total_volume, (1/3))) * volume_scalar) # Note this is the side of one length of the box

        ## NEED A CLAUSE TO ENSURE BOX length is bigger than the longest molecule ##

        packmol_input_filename = model.model_name + ".inp"
        packmol_input_filepath = os.path.join(model_dirs.packmol_inputs, packmol_input_filename)

        packmol_output_filename = model.model_name + ".pdb"
        packmol_output_filepath = os.path.join(model_dirs.packmol_systems, packmol_output_filename)
      
        lines = []
    
        if tolerance is None:
            tolerance = 2.0
        if filetype is None:
            filetype = "pdb"
        tolerance = "tolerance " + str(tolerance)
        ter = "add_amber_ter"
        output = "output " + packmol_output_filepath
        filetype = "filetype " + filetype
        lines.extend([tolerance, ter, output, filetype])

        def write_packmol_struct_block(model_dirs, lines, name, ratio, total_mols_in_model, size):#
         
            name = "structure " + model_dirs.molecules_dir + f"/{name}/{name}.pdb" 
            if ratio == False:
                number = "	" + "number " + str(1*molecule_scalar)
            else:
                number = "	" + "number " + str((int(total_mols_in_model*ratio))*molecule_scalar)    
            cube = "	inside cube 0. 0. 0. " + str(size+5) + "."
            end = "end structure"
            lines.extend([name, number, cube, end])
            return(lines) 

        # NOTE ALL RATIOS ARE NOW CALCULATED CORRECTLY and sum to 1.0 in the model object!! ##
        for i in range(len(model.molecules)):
            lines = write_packmol_struct_block(model_dirs, lines, model.molecules[i].name, model.molecule_ratios[i], model.min_mols_for_sim, size)
         
        #print(input_filepath)
        f = open(packmol_input_filepath, "w")
        for line in lines:
            f.write(line)
            f.write('\n')
        f.close()   

        packmol_command = str(model_dirs.packmol_path) + " < " + packmol_input_filepath
        print(packmol_command)
        try:
            result = subprocess.run(packmol_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
            # Command executed successfully
                print("Output:", result.stdout)
            else:
            # Command failed, print error message
                print("Error:", result.stderr)
        except Exception as e:
        # Exception occurred during subprocess execution
            print("Exception:", e)
    
        return()

    @staticmethod
    def extract_unique_rescodes(pdb_file):
        unique_residues = set()
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    residue_code = line[17:20].strip()
                    unique_residues.add(residue_code)
        return unique_residues

    @staticmethod
    def load_molecule_list(manager):
        molecule_dict = {}
        with open(manager.residue_code_csv, 'r') as file:
            for line in file:
                parts = line.strip().split(',')
                if len(parts) == 3:
                    molecule_dict[parts[2]] = parts[0]
        return molecule_dict

    @staticmethod
    def find_matching_molecules(residue_codes, molecule_dict):
        matching_molecules = []
        for code in residue_codes:
            if code in molecule_dict:
                matching_molecules.append(molecule_dict[code])
        return matching_molecules

    @staticmethod
    def generate_amber_params_from_packmol_bio_oil(manager, molecule_list, system_pdb_path):
        system_name = (system_pdb_path.split("/")[-1]).split(".")[0]

         #Write and exectute the intleap file. PROBLEMS LOADING SYSTEM PDB -not sure why though
        file_content = ""
        file_content = "source leaprc.gaff\n"

        for molecule in molecule_list:
            prepi_path = os.path.join(manager.molecules_dir, molecule, (molecule + ".prepi"))
            file_content += f"loadamberprep {prepi_path}\n"
            frcmod_path = os.path.join(manager.molecules_dir, molecule, (molecule + ".frcmod"))
            file_content += f"loadamberparams {frcmod_path}\n"

        file_content += f"system = loadpdb {system_pdb_path}\n"

        output_dir = os.path.join(manager.bio_oil_systems_dir, system_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
        prmtop_file = os.path.join(output_dir, (system_name + ".prmtop"))
        rst7_file = os.path.join(output_dir, (system_name + ".rst7"))


        #x, y, z = build.get_xyz_dists(system_pdb) # get x,y,z for periodic dists
        #x, y, z = x+2, y+2, z+2 # Add a buffer into it
        #file_content += f"setBox system vdw {x} {y} {z}\n" 

        file_content += f"setBox system centers\n" # this method of setting pbc didnt work. ACTUALLY DID WORK, just in the wrong place
        file_content += f"saveamberparm system {prmtop_file} {rst7_file}\n"

        file_content += "quit\n"

        intleap_path = os.path.join(output_dir, (system_name + ".intleap"))


        with open(intleap_path, 'w') as file:
            file.write(file_content)
            
        leap_command = "tleap -f " + intleap_path
        #print("The command that would be run in the shell is: ")
        #print(leap_command)
        print(intleap_path)
        try:
            result = subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
                # Command executed successfully
                print("Output:", result.stdout)
            else:
                 # Command failed, print error message
                print("Error:", result.stderr)
        except Exception as e:
                     # Exception occurred during subprocess execution
                    print("Exception:", e)
    