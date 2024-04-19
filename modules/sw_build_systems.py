# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:00:30 2024

@author: danie
"""
import csv
import itertools
import os
import subprocess
import numpy as np
import shutil
from openbabel import openbabel, pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolFromPDBFile

class BuildSystems():
    def __init__(self):
        pass 
    
    def SmilesToPDB(self, smiles_string, output_file):
        """
        Converts a SMILES string to a PDB file.

        Args:
            smiles_string (str): The SMILES string of the molecule.
            output_file (str): The name of the output PDB file.

        Returns:
            None. Writes the 3D structure of the molecule to a PDB file.
            
        Note: 
            This function is utilised by the SmilesToPDB_GenerateRescode function and carries out
            the same functionality but additionally generates a residue code for the pdb file generated.
            These generated residue codes are stored in a database.
        """
        # Create a molecule from the SMILES string
        molecule = pybel.readstring("smi", smiles_string)

        # Generate the 3D structure
        molecule.make3D()

        # Convert the molecule to PDB format
        pdb_string = molecule.write("pdb")

        # Write the PDB string to a file
        with open(output_file, "w") as f:
            f.write(pdb_string)

    def SmilesToPDB_GenerateRescode(self, smiles, name, directories):
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
        residue_codes = self.load_residue_codes(directories.residue_code_csv)
        # Check if the name or smiles is already in the database
        existing_entry = self.find_existing_entry(residue_codes, name, smiles)
        if existing_entry:
            # Use existing residue code
            residue_code = existing_entry[2]  # Assuming code is the third column
        else:
            # Generate a unique 3-letter residue code excluding forbidden codes
            residue_code = self.generate_unique_residue_code(residue_codes, forbidden_codes)
            # Update the CSV file with the new entry
            self.update_residue_codes_csv(name, smiles, residue_code, directories.residue_code_csv)
        # Replace default "UNL" codes in the PDB content - NOTE, "UNL" is also a forbidden code as it is the default code.
        pdb_filepath = os.path.join(directories.pdb_file_dir, (name + ".pdb"))
        self.SmilesToPDB(smiles, pdb_filepath)    
        with open(pdb_filepath, "r") as pdb_file:
            lines = pdb_file.readlines()     
            for i, line in enumerate(lines):
                lines[i] = line.replace(" UNL", f" {residue_code}")   
        with open(pdb_filepath, "w") as pdb_file:              
            pdb_file.writelines(lines)
        return()

    def load_residue_codes(self, residue_code_csv):
        """
        Loads existing residue codes from a CSV file. Used by SmilesToPDB_GenerateRescode function.

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

    def generate_unique_residue_code(self, residue_codes, forbidden_codes=None):
        """
        Generates a unique 3-letter residue code not present in the database and not in the forbidden codes list.
        Used by SmilesToPDB_GenerateRescode function.

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

    def update_residue_codes_csv(self, name, smiles, residue_code, residue_code_csv):
        """
        Updates a CSV file with a new entry if the entry does not already exist.
        Used by SmilesToPDB_GenerateRescode function.

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
        residue_codes = self.load_residue_codes(residue_code_csv)
        existing_entry = self.find_existing_entry(residue_codes, name, smiles)

        if not existing_entry:
            with open(residue_code_csv, "a", newline="") as csv_file:
                writer = csv.writer(csv_file)

                # Add a new entry to the CSV file
                writer.writerow([name, smiles, residue_code])

    def find_existing_entry(self, residue_codes, name, smiles):
        """
        Finds an existing entry in a list of residue codes based on molecule name or SMILES representation.
        Used by SmilesToPDB_GenerateRescode function.

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

    def max_pairwise_distance(self, mol):
        """
        Calculates the maximum pairwise distance between atoms in a molecule.

        Parameters:
            - mol (rdkit.Chem.Mol): RDKit molecule object.

        Returns:
            - float: Maximum pairwise distance between 2 atoms in the molecule.
        """
        conformer = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        # Extract atom coordinates
        atom_positions = np.zeros((num_atoms, 3))
        for i in range(num_atoms):
            pos = conformer.GetAtomPosition(i)
            atom_positions[i] = (pos.x, pos.y, pos.z)
        # Calculate pairwise distances
        distances = np.linalg.norm(atom_positions[:, np.newaxis, :] - atom_positions, axis=2)
        # Exclude self-distances and get the maximum distance
        np.fill_diagonal(distances, 0)
        max_distance = np.max(distances)
        return max_distance

class BuildAmberSystems(BuildSystems):
    
    error_param = "Molecule not parametrized. Please parametrize pdb_file."
    
    def __init__(self):
        pass
    
    def mod_pdb_file(self, pdb_file): # pdb_file is the entire path to the pdb_file - see jupyter notebook
        """
        Modify a PDB file by removing duplicate bond information.
         Note: parametrization does not work if double bond connections are specified.

        Parameters:
        - pdb_file (str): The full path to the PDB file.

        Returns:
        None

        This function reads a PDB file, removes duplicate bond information in lines starting with "CONECT",
        and writes the modified content back to the same file.

        Note: The modification involves removing duplicate atom connection information while maintaining the original order.
        """
        file_path = pdb_file  
        # Open the file and read its lines
        with open(file_path, 'r') as file:
            lines = file.readlines()
        # Process each line
        modified_lines = []
        for line in lines:
            # Process lines starting with "CONECT"
            if line.startswith("CONECT"):
                # Split the line into numbers
                numbers = line.split()[1:]
                # Remove duplicates while maintaining order
                unique_numbers = sorted(set(numbers), key=numbers.index)
                # Join the modified numbers and reconstruct the line
                modified_line = f"{line.split()[0]:<8}" + " ".join(unique_numbers) + "\n"
                modified_lines.append(modified_line)
            else:
                modified_lines.append(line)
        # Write the modified lines back to the file
        with open(file_path, 'w') as file:
            file.writelines(modified_lines)
    
    def parametrize_mol(self, directories, molecule_name):
        # Create a new directory for param files for the molecule and copy pdb there
        pdb_filepath = os.path.join(directories.pdb_file_dir, (molecule_name + ".pdb"))
        self.mod_pdb_file(pdb_filepath)
        param_mol_dir = os.path.join(directories.molecules_dir, molecule_name)
        if not os.path.exists(param_mol_dir):
            os.makedirs(param_mol_dir, exist_ok=True)
        shutil.copy2(pdb_filepath, param_mol_dir)
        os.remove(pdb_filepath)
        
        # Specify paths for tleap
        pdb_filepath = os.path.join(param_mol_dir, (molecule_name + ".pdb"))
        mol2_filepath = os.path.join(param_mol_dir, (molecule_name + ".mol2"))
        antechamber_command = "antechamber -i " + pdb_filepath + " -fi pdb -o " + mol2_filepath + " -fo mol2 -c bcc -s 2"       
        #print("The antechamber command that will run is:")
        #print(antechamber_command)
        subprocess.run(antechamber_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
        frcmod_filepath = os.path.join(param_mol_dir, (molecule_name + ".frcmod"))
        parmchk_command = "parmchk2 -i " + mol2_filepath + " -f mol2 -o " + frcmod_filepath
        #print("The parmchk command that will run is:")
        #print(parmchk_command)      
        subprocess.run(parmchk_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        lib_filepath = os.path.join(param_mol_dir, (molecule_name + ".lib"))
        intleap_filepath = os.path.join(param_mol_dir, (molecule_name + ".intleap"))
        
        file_content = f"""source leaprc.protein.ff14SB
        source leaprc.gaff
        source leaprc.water.fb3
        {molecule_name} = loadmol2 {mol2_filepath}
        saveoff {molecule_name} {lib_filepath}
        quit
        """
        with open(intleap_filepath, 'w') as file:
            file.write(file_content)
        leap_command = "tleap -f " + intleap_filepath
        #print("The leap command that will run is:")
        #print(leap_command)
        subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
   
    def is_mol_parametrized(self, directories, molecule_name):
        param_mol_dir = os.path.join(directories.molecules_dir, molecule_name)
        molecule_dir = os.path.join(directories.pdb_file_dir, molecule_name)
        pdb_filepath = os.path.join(param_mol_dir, (molecule_name + ".pdb"))
        mol2_filepath = os.path.join(param_mol_dir, (molecule_name + ".mol2"))
        files_exist = os.path.exists(pdb_filepath) and os.path.exists(mol2_filepath)
        return(files_exist) # This will be true or false
    
    def gen_3_3_array(self, directories, molecule_name):
        """
        Generate a 3x3 array of molecules around a central molecule and solvate them in water.

        Parameters:
            - directories (object): Object containing directory paths.
            - molecule_name (str): Name of the central molecule.

        Returns:
            None
            
        Note:
            Nothing is returned - but a new set of files for an amber simulation array is generated
        """
        if self.is_mol_parametrized(directories, molecule_name) == True:
            pass
        if self.is_mol_parametrized(directories, molecule_name) == False:
            print(self.error_param)
            return()
        
        file_subtype = "_3_3_array"
        
        pdb_filepath = os.path.join(directories.molecules_dir, molecule_name, (molecule_name + ".pdb"))
        mol2_filepath = os.path.join(directories.molecules_dir, molecule_name, (molecule_name + ".mol2"))
        output_dir = os.path.join(directories.systems_dir, (molecule_name + file_subtype))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        mol = MolFromPDBFile(pdb_filepath)
        max_dist = self.max_pairwise_distance(mol)
        translate_distance = float((int(max_dist)+1))
        box_dist = float((int(max_dist)+1)*3)
        
        molecule_name_1 = molecule_name + "_1"
        molecule_name_2 = molecule_name + "_2"
        molecule_name_3 = molecule_name + "_3"
        molecule_name_4 = molecule_name + "_4"
        molecule_name_5 = molecule_name + "_5"
        molecule_name_6 = molecule_name + "_6"
        molecule_name_7 = molecule_name + "_7"
        molecule_name_8 = molecule_name + "_8"
        molecule_name_9 = molecule_name + "_9"

        translate_line_1 = "{0.0 0.0 0.0}"
        translate_line_2 = "{0.0 0.0 " + str(translate_distance) + "}"
        translate_line_3 = "{0.0 0.0 " + str(-translate_distance) + "}"

        translate_line_4 = "{0.0 " + str(translate_distance) + " " + str(translate_distance) + "}"
        translate_line_5 = "{0.0 " + str(translate_distance) + " " + str(-translate_distance) + "}"
        translate_line_6 = "{0.0 " + str(translate_distance) + " 0.0}"

        translate_line_7 = "{0.0 " + str(-translate_distance) + " " + str(translate_distance) + "}"
        translate_line_8 = "{0.0 " + str(-translate_distance) + " " + str(-translate_distance) + "}"
        translate_line_9 = "{0.0 " + str(-translate_distance) + " 0.0}"

        combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + " " + molecule_name_5 + " " + molecule_name_6 + " " + molecule_name_7 + " " + molecule_name_8 + " " + molecule_name_9 + "}" 

        intleap_path = os.path.join(output_dir, (molecule_name + file_subtype + ".intleap"))
        prmtop_filepath =  os.path.join(output_dir, molecule_name + file_subtype + ".prmtop")
        rst_filepath = os.path.join(output_dir, molecule_name + file_subtype + ".rst7")
        three_three_array_pdb_filepath = os.path.join(output_dir, molecule_name + file_subtype + ".pdb")

        file_content = f"""source leaprc.protein.ff14SB
        source leaprc.gaff
        source leaprc.water.fb3
        {molecule_name_1} = loadMol2 {mol2_filepath}
        {molecule_name_2} = loadMol2 {mol2_filepath}
        {molecule_name_3} = loadMol2 {mol2_filepath}
        {molecule_name_4} = loadMol2 {mol2_filepath}
        {molecule_name_5} = loadMol2 {mol2_filepath}
        {molecule_name_6} = loadMol2 {mol2_filepath}
        {molecule_name_7} = loadMol2 {mol2_filepath}
        {molecule_name_8} = loadMol2 {mol2_filepath}
        {molecule_name_9} = loadMol2 {mol2_filepath}

        translate {molecule_name_1} {translate_line_1}
        translate {molecule_name_2} {translate_line_2}
        translate {molecule_name_3} {translate_line_3}
        translate {molecule_name_4} {translate_line_4}
        translate {molecule_name_5} {translate_line_5}
        translate {molecule_name_6} {translate_line_6}
        translate {molecule_name_7} {translate_line_7}
        translate {molecule_name_8} {translate_line_8}
        translate {molecule_name_9} {translate_line_9}
         
        system = combine {combine_line}

        solvateBox system TIP3PBOX {box_dist}

        saveamberparm system {prmtop_filepath} {rst_filepath}
        savepdb system {three_three_array_pdb_filepath}
        quit
        """
        with open(intleap_path, 'w') as file:
            file.write(file_content)
            
        leap_command = "tleap -f " + intleap_path
        #print("The command that would be run in the shell is: ")
        #print(leap_command)
        subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    def gen_2_2_array(self, directories, molecule_name):
        """
        Generate a 2x2 array of molecules and solvate them in water.

        Parameters:
            - directories (object): Object containing directory paths.
            - molecule_name (str): Name of the central molecule.

        Returns:
            None
            
        Note:
            Nothing is returned - but a new set of files for an amber simulation array is generated
        """
        if self.is_mol_parametrized(directories, molecule_name) == True:
            pass
        if self.is_mol_parametrized(directories, molecule_name) == False:
            print(self.error_param)
            return()
        
        file_subtype = "_2_2_array"
        
        pdb_filepath = os.path.join(directories.molecules_dir, molecule_name, (molecule_name + ".pdb"))
        mol2_filepath = os.path.join(directories.molecules_dir, molecule_name, (molecule_name + ".mol2"))
        output_dir = os.path.join(directories.systems_dir, (molecule_name + file_subtype))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        mol = MolFromPDBFile(pdb_filepath)
        max_dist = self.max_pairwise_distance(mol)

        # We will translate the molecules from the centre by their max pariwise distance.
        # This will give an overall distance of (2*max_pairwise_distance) between each molecule
        translate_dist = float((int(max_dist)+1))

        # The box_distance is (4*max_pairwise_distance) as we have 4 molecules
        box_dist = float((int(max_dist)+1)*4) # the +1 is here as the "int" function converts a decimal to an integer but always rounds down
    
        # Need to load four instances of the same molecule
        molecule_name_1 = molecule_name + "_1"
        molecule_name_2 = molecule_name + "_2"
        molecule_name_3 = molecule_name + "_3"
        molecule_name_4 = molecule_name + "_4"

        # Need four individual translate lines
        translate_line_1 = "{0.0 " + str(translate_dist) + " " + str(translate_dist) + "}"
        translate_line_2 = "{0.0 " + str(-translate_dist) + " " + str(translate_dist) + "}"
        translate_line_3 = "{0.0 " + str(translate_dist) + " " + str(-translate_dist) + "}"
        translate_line_4 = "{0.0 " + str(-translate_dist) + " " + str(-translate_dist) + "}"

        # Need to create a system of our four instances of the molecules
        combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + "}" 

        intleap_path = os.path.join(output_dir, (molecule_name + file_subtype + ".intleap"))
        prmtop_filepath =  os.path.join(output_dir, molecule_name + file_subtype + ".prmtop")
        rst_filepath = os.path.join(output_dir, molecule_name + file_subtype + ".rst7")
        two_two_array_pdb_filepath = os.path.join(output_dir, molecule_name + file_subtype + ".pdb")
        
        # File content for the intleap file with all of the variables specified above
        file_content = f"""source leaprc.protein.ff14SB
        source leaprc.gaff
        source leaprc.water.fb3
        {molecule_name_1} = loadMol2 {mol2_filepath}
        {molecule_name_2} = loadMol2 {mol2_filepath}
        {molecule_name_3} = loadMol2 {mol2_filepath}
        {molecule_name_4} = loadMol2 {mol2_filepath}

        translate {molecule_name_1} {translate_line_1}
        translate {molecule_name_2} {translate_line_2}
        translate {molecule_name_3} {translate_line_3}
        translate {molecule_name_4} {translate_line_4}
         
        system = combine {combine_line}

        saveamberparm system {prmtop_filepath} {rst_filepath}
        savepdb system {two_two_array_pdb_filepath}
        quit
        """
        # Write content to the intleap filepath
        with open(intleap_path, 'w') as file:
            file.write(file_content)
    
        leap_command = "tleap -f " + intleap_path
        #print("The command that would be run in the shell is: ")
        #print(leap_command)
        subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)