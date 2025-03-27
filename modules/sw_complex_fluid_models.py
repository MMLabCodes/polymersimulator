from modules.sw_orca import *
from modules.sw_basic_functions import *
import os as os
import subprocess
import csv
import math
import pandas as pd

class complex_fluid_model:
    # This class is the output from the functions that make models and allows for a standardised version of them
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
    # This  class contains all the functions to generate complex_fluid_models as described in the paper
    # The inputs are orca molecules
    # The outputs are instances of the 'complex_fluid_model' class - all standardised so it is easier to work with
    
    def __init__(self):
        pass

    @staticmethod
    def get_weighted_average(complex_fluid_model, class_attribute): # class attribute and model type must be passed as a "string"
        peak_area = []
        weighted_average = 0

        # Not sure what this function is, there is already a function developed to calculat eheteroatom content of things    
        if "Heteroatom_content" in class_attribute:
            # Everything after the else for working with floats or integers that are attributes of the molecule class
            # Below we work out the weighted average of the heteroatoms content
            def get_heteroatom_content(smiles):
                heteroatom_wt = 0
                if class_attribute == "Heteroatom_content_O":
                    het = smiles.count("O")
                    het = het + smiles.count("o")
                    heteroatom_wt = het * 16
                    return(heteroatom_wt)
                if class_attribute == "Heteroatom_content_N":
                    het = smiles.count("N")
                    het = het + smiles.count("n")
                    heteroatom_wt = het * 14
                    return(heteroatom_wt)
                if class_attribute == "Heteroatom_content_S":
                    het = smiles.count("S")
                    het = het + smiles.count("s")
                    heteroatom_wt = het * 32
                    return(heteroatom_wt)
             
                return(return_list)
        
            weighted_het_weight = 0
            weighted_het_weight = sum(ratio * get_heteroatom_content(molecule.smiles) for ratio, molecule in zip(complex_fluid_model.molecule_ratios, complex_fluid_model.molecules))
            return(weighted_het_weight)    
        else:
            weighted_average = sum(ratio * float(getattr(molecule, class_attribute)) for ratio, molecule in zip(complex_fluid_model.molecule_ratios, complex_fluid_model.molecules))
                                         
            return(weighted_average)
            
    @staticmethod
    def calculate_heteroatom_percentages_sing_mol(molecule):
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
        def get_mixture(molecules): # model is a a list of python objects   
            mixture = {obj.smiles: float(obj.peak_area) for obj in molecules}
            return(mixture)
            
        mixture = get_mixture(molecules) 
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
        Check if the given molecule is an instance of the orca_molecule class.

        Parameters:
        -----------
        molecule : object
        The object to check.

        Returns:
        --------
        bool
        True if the object is an instance of orca_molecule, False otherwise.
        """
        return isinstance(molecule, orca_molecule)

    @staticmethod
    def all_model(model_name, orca_molecules):
        orca_mol_bool = [isinstance(molecule, orca_molecule) for molecule in orca_molecules]
        if "False" in orca_mol_bool:
            print("Please use orca molecule class instances for each molecule in list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return(None)

        ratios = [float(mol.peak_area)/100 for mol in orca_molecules]
        all_model = complex_fluid_model(orca_molecules, None, ratios, "ALL_model", model_name)
        
        return(all_model)

    @staticmethod
    def fixed_threshold_model(model_name, orca_molecules, selection_threshold):
        orca_mol_bool = [isinstance(molecule, orca_molecule) for molecule in orca_molecules]
        if "False" in orca_mol_bool:
            print("Please use orca molecule class instances for each molecule in list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return(None)

        fixed_threshold_model = [] 

        for molecule in orca_molecules:
            if float(molecule.peak_area) >= selection_threshold:
                fixed_threshold_model.append(molecule)

        tot_peak_area = sum([float(mol.peak_area) for mol in fixed_threshold_model])
        ratios = [(float(mol.peak_area)/tot_peak_area) for mol in fixed_threshold_model]
        
        fixed_threshold_model = complex_fluid_model(fixed_threshold_model, None, ratios, "FT_model", model_name)
        
        return(fixed_threshold_model)

    @staticmethod
    def proportional_threshold_model(model_name, orca_molecules):
        orca_mol_bool = [isinstance(molecule, orca_molecule) for molecule in orca_molecules]
        if "False" in orca_mol_bool:
            print("Please use orca molecule class instances for each molecule in list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return(None)

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
        # used in both AG and SG models
        tot_peak_area = 0
        for molecule in orca_molecules:
            tot_peak_area = tot_peak_area + float(molecule.peak_area)
        return(tot_peak_area)
        
    @staticmethod
    def abundancy_grouped_model(model_name, orca_molecules):
        orca_mol_bool = [isinstance(molecule, orca_molecule) for molecule in orca_molecules]
        if "False" in orca_mol_bool:
            print("Please use orca molecule class instances for each molecule in list passed to this function.")
            print("Refer to the sw_orca module for more information on setting this up.")
            return(None)

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
    def scored_grouped_model(model_name, orca_molecules):
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
    def rank_models(df):
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
    def min_mols_4_simulation(model):
        num_compounds = len(model.molecules)
        min_peak = min(model.molecule_ratios)
        min_mols = 0
        for i in range(num_compounds):
            min_mols += (float(model.molecule_ratios[i]) / float(min_peak)) * (len(model.molecules) + 1)
    
        return int(min_mols)

    @staticmethod
    def min_atoms_4_simulation(model):
        min_peak = min(model.molecule_ratios)
        num_atoms = 0
        for i in range(len(model.molecules)):
            mol = Chem.MolFromSmiles(model.molecules[i].smiles)
            mol = Chem.AddHs(mol)
            num_atoms = num_atoms + (mol.GetNumAtoms()*(float(model.molecule_ratios[i]) / float(min_peak))*(len(model.molecules) + 1))
            
        return int(num_atoms)
        
    @staticmethod
    def min_vol_4_simulation(model):
        tot_volume = 0      
        for i in range(len(model.molecules)):
            volume_of_mol = model.molecules[i].volume*(model.molecule_ratios[i]*model.min_mols_for_sim)
            tot_volume = tot_volume + volume_of_mol
         
        return tot_volume

    @staticmethod
    def model_output_block(all_model, model, ranked_data):
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
    