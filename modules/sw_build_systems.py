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

    packmol_path = "/home/s.983045/bin/packmol/packmol-20.14.0/packmol"
    
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

    def PolymerUnits_GenerateRescode(self, directories, name):
        # Name should be of a trimer
        if "trimer" not in name:
            print("Polymeric unit generation requires trimers. Please consult the build systems guide for information on how to do this")
            return()
    
        forbidden_codes = ["AAA", "BBB", "CCC", "UNL", "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "SEC", "TRP", "TYR", "VAL"]
   
        # Load existing residue codes from the CSV file
        residue_codes = self.load_residue_codes(directories.residue_code_csv)
    
        # Check if the name or smiles is already in the database
        existing_entry = self.find_existing_entry(residue_codes, name)
        if existing_entry:
            # Use existing residue code
            residue_code = existing_entry[2]  # Assuming code is the third column  
            residue_smiles = existing_entry[1] # Assuming SMILES is the second column

        if not existing_entry:
            print("Please parameterize the initial trimer unit and generate a residue code for it using it's SMILES string")
            return()

        # New codes and names for residue code csv are specified here.
        head_code = "h" + residue_code[1:]
        mainchain_code = "m" + residue_code[1:]
        tail_code = "t" + residue_code[1:]

        head_name = "head_" + name
        mainchain_name = "mainchain_" + name
        tail_name = "tail_" + name

        # A tag is also placed in front of the SMILES. Note, this should not be required, but the residue code file will not update unless a new SMILES AND molecule name are passed to it. In this case, the SMILES for each unit are simply the same SMILES as the trimer but with a h,m,t tag in front of the SMILES.
        head_smiles = "h" + residue_smiles
        mainchain_smiles = "m" + residue_smiles
        tail_smiles = "t" + residue_smiles
    
        self.update_residue_codes_csv(head_name, head_smiles, head_code, directories.residue_code_csv)
        self.update_residue_codes_csv(mainchain_name, mainchain_smiles, mainchain_code, directories.residue_code_csv)
        self.update_residue_codes_csv(tail_name, tail_smiles, tail_code, directories.residue_code_csv)
        print("Head code assigned: ", head_code)
        print("Mainchain code assigned: ", mainchain_code)
        print("Tail code assigned: ", tail_code)
    
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

    def find_existing_entry(self, residue_codes, name, smiles=None):
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
        if smiles == None:
            for entry in residue_codes:
                if len(entry) >= 2 and entry[0] == name:
                    return entry
        else:
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

        Notes:
            This function will no longer be updated.
        """
        warnings.warn(
            "The max_pairwise_distance method is deprecated and will be removed in a future version. "
            "Please use the get_xyz_dists method instead.",
            DeprecationWarning
        )
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

    def run_packmol(self, directories, input_file_name):
        """
        Run Packmol with the specified input file and directories.

        Args:
            directories: An object with methods to load Packmol file paths and access system directories.
            input_file_name: The name of the Packmol input file.

        The 'directories' object should have the following methods:
            - load_pckml_filepath(input_file_name): Returns the full path to the Packmol input file.
            - systems_dir: A directory where the Packmol output will be stored.

        Example:
            from modules.sw_directories import *
            directories = PolymerSimulatorDirs("path/to/project/directory")

            from modules.sw_build_systems import *
            build = BuildSystems()
            build.run_packmol(directories, 'input_file_name')
        """
        if os.path.exists(self.packmol_path):
            print(f"Packmol executable exists at '{self.packmol_path}'.")
        else:
            print(f"Packmol executable not found at '{self.packmol_path}'.")
            print("Please update the class variable for the packmol path using 'update_packmol_path(new_packmol_path)'")
            print("")
            print("Call packmol_help() for more information")
            return()
        packmol_filepath = directories.load_pckml_filepath(input_file_name)
        if os.path.exists(packmol_filepath):
            print(f"Packmol input file exists at '{packmol_filepath}'.")
        else:
            print("Packmol input file not found.")
            print("")
            return()
        system_dir = os.path.join(directories.systems_dir, input_file_name)
        packmol_command = self.packmol_path + " < " + packmol_filepath
        result = subprocess.run(packmol_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return()

    @classmethod
    def update_packmol_path(cls, new_packmol_path):
        """
        Update the path to the Packmol executable.

        Args:
            new_packmol_path: The new path to the Packmol executable.

        Example:
            from modules.sw_build_systems import *
            build = BuildSystems()
            build.update_packmol_path('/new/path/to/packmol/ececutable')
        """
        if os.path.exists(new_packmol_path):
            self.packmol_path = new_packmol_path
            print(f"Packmol path updated to '{self.packmol_path}'.")
        else:
            print(f"The provided path '{new_packmol_path}' does not exist. Please provide a valid path.")

    @staticmethod
    def get_xyz_dists(input_file=None):
        """
        Calculates the maximum distance between the largest and smallest xyz coordinates.

        Parameters:
            - input_file: Currently supported formats are '.pdb' and '.xyz' files.

        Returns:
            - tuple: Maximum distances along x, y, and z axes between 2 atoms.

        Example:
            from modules.sw_build_systems import *
            build = BuildSystems()
            pdb_file = path/to/your/pdb/file
            x, y, z = build.get_xyz_dists(pdb_file)
        """
        if input_file == None:
            print("Please provide a valid input file.")
            print("Supported formats are: '.pdb' and '.xyz'")
            print("")
            #self.get_xyz_dists_help()
            return()
        
        if ".pdb" in input_file:
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat("pdb")
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, input_file)
            # Extract atom coordinates
            coords = []
            for atom in openbabel.OBMolAtomIter(mol):
                coords.append((atom.GetX(), atom.GetY(), atom.GetZ()))

        if ".xyz" in input_file:
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat("xyz")
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, file_path)
            coords = []
            for atom in openbabel.OBMolAtomIter(mol):
                coords.append((atom.GetX(), atom.GetY(), atom.GetZ()))
    
        # Convert to numpy array for easy manipulation
        coords = np.array(coords)
        
        # Calculate max and min for each axis
        max_x, max_y, max_z = np.max(coords, axis=0)
        min_x, min_y, min_z = np.min(coords, axis=0)

        # Calculate the distances
        dist_x = max_x - min_x
        dist_y = max_y - min_y
        dist_z = max_z - min_z

        return(dist_x, dist_y, dist_z)

        @classmethod
        def get_xyz_dists_help(cls):
            """Display help information for the get_xyz_dists method."""
            print(cls.get_xyz_dists.__doc__)
            print("This method calculates the maximum distance between the largest and smallest xyz coordinates")
            print("from a PDB or XYZ file.")

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
    
    def parametrize_mol(self, directories=None, molecule_name=None):
        if directories == None or molecule_name == None:
            print("Please provide 2 arguments as follows: parametrize_mol(directories, molecule_name)")
            print("Directories: A python object generated with the PolymerSimulatorDirs(filepath) method imported from sw_directories")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            return(None)
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

    def gen_ac_file(self, directories=None, molecule_name=None):
        if directories == None or molecule_name == None:
            print("Please provide 2 arguments as follows: gen_ac_file(directories, molecule_name)")
            print("Directories: A python object generated with the PolymerSimulatorDirs(filepath) method imported from sw_directories")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            return(None)
        mol2_name = os.path.join(directories.molecules_dir, molecule_name, molecule_name) + ".mol2"
        ac_name = os.path.join(directories.molecules_dir, molecule_name, molecule_name) + ".ac"
        antechamber_command = "antechamber -fi mol2 -fo ac -i " + mol2_name + " -o " + ac_name + " -c bcc -s 2"
        subprocess.run(antechamber_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return(None)

    def gen_prepin_files(self, directories, molecule_name):
        if directories == None or molecule_name == None:
            print("Please provide 2 arguments as follows: gen_prepin_files(directories, molecule_name)")
            print("Directories: A python object generated with the PolymerSimulatorDirs(filepath) method imported from sw_directories")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            return(None)
        
        # Find the path of the molecule directory (where all the required files are required)
        molecule_dir = os.path.join(directories.molecules_dir, molecule_name)
        if os.path.exists(molecule_dir):
            pass
        else:
            print("Please parameterize molecule.")
            return()  
    
        # Define file names required for prepgen
        # Could have an additional function that generates these files
        head_file = "head_" + molecule_name + ".txt"
        head_prepin_file = os.path.join(molecule_dir, "head_" + molecule_name + ".prepi")
        mainchain_file = "mainchain_" + molecule_name + ".txt"
        mainchain_prepin_file = os.path.join(molecule_dir, "mainchain_" + molecule_name + ".prepi")
        tail_file = "tail_" + molecule_name + ".txt"
        tail_prepin_file = os.path.join(molecule_dir, "tail_" + molecule_name + ".prepi")
        ac_name = os.path.join(molecule_dir, molecule_name + ".ac")
    
        # Check if the required files exist
        segment_files = [head_file, mainchain_file, tail_file]
        for file_path in segment_files:
            if os.path.exists(os.path.join(molecule_dir, file_path)):
                pass
            else:
                print("Please prepare files defining each segment of the trimer. Head, mainchain and tail. Call build_systems.prepin_help() for more information.")
                return() 

        # Retrieve the residue codes from the database
        head_rescode, mainchain_rescode, tail_rescode = directories.retrieve_polymeric_rescodes("3HB_trimer")
        print(head_rescode, mainchain_rescode, tail_rescode)
        if head_rescode == None or mainchain_rescode == None or tail_rescode == None:
            print("Residue codes for polymeric units not generated. Please generate them with 'build.PolymerUnits_GenerateRescode(directories, 'ethane')'")
            return()

        # Define the prepgen commands with the collated information
        head_command = "prepgen -i " + ac_name + " -o " + head_prepin_file + " -f prepi -m " + head_file + " -rn " + head_rescode + " -rf " + head_rescode + ".res"
        mainchain_command = "prepgen -i " + ac_name + " -o " + mainchain_prepin_file + " -f prepi -m " + mainchain_file + " -rn " + mainchain_rescode + " -rf " + mainchain_rescode + ".res"
        tail_command = "prepgen -i " + ac_name + " -o " + tail_prepin_file + " -f prepi -m " + tail_file + " -rn " + tail_rescode + " -rf " + tail_rescode + ".res"
        commands = [head_command, mainchain_command, tail_command]
        os.chdir(molecule_dir)
        for command in commands:
            try:
                result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if result.returncode == 0:
                    # Command executed successfully
                    print("Output:", result.stdout)
                else:
                    # Command failed, print error message
                    print("Error:", result.stderr)
            except Exception as e:
                # Exception occurred during subprocess execution
                print("Exception:", e)
        os.chdir(directories.main_dir)
        return(None)
    
    def is_mol_parametrized(self, directories, molecule_name):
        param_mol_dir = os.path.join(directories.molecules_dir, molecule_name)
        molecule_dir = os.path.join(directories.pdb_file_dir, molecule_name)
        pdb_filepath = os.path.join(param_mol_dir, (molecule_name + ".pdb"))
        mol2_filepath = os.path.join(param_mol_dir, (molecule_name + ".mol2"))
        files_exist = os.path.exists(pdb_filepath) and os.path.exists(mol2_filepath)
        return(files_exist) # This will be true or false

    def is_poly_prepped(self, directories, molecule_name):
        poly_prepped_dir = os.path.join(directories.molecules_dir, molecule_name)
        head_prepi_filepath = os.path.join(poly_prepped_dir, ("head_" + molecule_name + ".prepi"))
        mainchain_prepi_filepath = os.path.join(poly_prepped_dir, ("mainchain_" + molecule_name + ".prepi"))
        tail_prepi_filepath = os.path.join(poly_prepped_dir, ("tail_" + molecule_name + ".prepi"))
        files_exist = os.path.exists(head_prepi_filepath) and os.path.exists(mainchain_prepi_filepath) and os.path.exists(tail_prepi_filepath)
        return(files_exist) # This will be true or false

    def gen_polymer_pdb(self, dirs, molecule_name, number_of_units):
        """
        Generates a polymer PDB file using `tleap` based on the specified molecule and number of units.

        Args:
            dirs (object): An object containing directory paths for the molecules and systems. Must have attributes
                       `molecules_dir` and `systems_dir`.
            molecule_name (str): The base name of the molecule to be used for generating the polymer.
            number_of_units (int): The number of units in the polymer.

        NOTE: The molecule name should be something like "3HB_trimer" and '.prepi' files should be available for that molecule.
            (these are required for using tleap to make polymers)

        Description:
            This function performs the following steps:
            1. Changes the working directory to the specified molecule's directory.
            2. Constructs file paths for the required input files (`head`, `mainchain`, and `tail` prepi files) and the
               output directory.
            3. Creates the output directory if it does not already exist.
            4. Constructs the polymer name and the paths for the output files (`prmtop`, `rst7`, and `pdb`).
            5. Retrieves residue codes for the polymeric units from the `dirs` object.
            6. Constructs the polymer sequence command for `tleap`.
            7. Creates the `tleap` input script (`.intleap` file) with the appropriate commands to generate the polymer.
            8. Executes the `tleap` command to generate the polymer, capturing and printing the output or errors.

        Raises:
            Exception: If there is an error changing the directory or executing the `tleap` command.

        Example Usage:
            dirs = DirectoryPaths('path/to/main/project/directory')
            gen_polymer_pdb(dirs, "3HB_trimer", 10)
        """
        molecule_dir = os.path.join(dirs.molecules_dir, molecule_name)
        try:
            os.chdir(molecule_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
        
        file_subtype = "_" + str(number_of_units) + "_polymer"
        head_prepi_filepath = "head_" + molecule_name + ".prepi"
        tail_prepi_filepath = "tail_" + molecule_name + ".prepi"
        mainchain_prepi_filepath = "mainchain_" + molecule_name + ".prepi"

        output_dir = os.path.join(dirs.systems_dir, (molecule_name.split("_")[0] + file_subtype))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        polymer_name = molecule_name.split("_")[0] + "_" + str(number_of_units) + "_polymer"

        intleap_path = polymer_name + ".intleap"
        prmtop_filepath = os.path.join(output_dir, polymer_name + ".prmtop")
        rst_filepath = os.path.join(output_dir, polymer_name + ".rst7")
        pdb_filepath = os.path.join(output_dir, polymer_name + ".pdb")

        head_rescode, mainchain_rescode, tail_rescode = dirs.retrieve_polymeric_rescodes(molecule_name)

        polymer_code = " ".join([head_rescode] + [mainchain_rescode] * (number_of_units - 2) + [tail_rescode])
        polymer_command = "{" + polymer_code + "}"

        file_content = f"""source leaprc.gaff
             source leaprc.water.fb3
             source leaprc.protein.ff14SB

             loadamberprep {head_prepi_filepath}
             loadamberprep {mainchain_prepi_filepath}
             loadamberprep {tail_prepi_filepath}

             list

             polymer = sequence {polymer_command}
             saveamberparm polymer {prmtop_filepath} {rst_filepath}
             savepdb polymer {pdb_filepath}
             quit
             """
        with open(intleap_path, 'w') as file:
             file.write(file_content)

        leap_command = "tleap -f " + intleap_path

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
        
        try:
            os.chdir(dirs.main_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
        
    
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
        unsolved_three_three_array_pdb_filepath = os.path.join(output_dir, "unsovled_" + molecule_name + file_subtype + ".pdb")

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

        savePDB system {unsolved_three_three_array_pdb_filepath}
        
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

    def build_3_3_polymer_array(self, directories=None, molecule_name=None, number_of_units=None):
         if directories == None or molecule_name == None or number_of_units == None:
            print("Please provide 3 arguments as follows: build_3_3_polymer_array(directories, molecule_name, number_of_units)")
            print("Directories: A python object generated with the PolymerSimulatorDirs(filepath) method imported from sw_directories")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            print("Number of units: An integer denoting polymer length, i.e. 10. NOTE: number_of_units >= 3")
            return(None)
        
         if self.is_poly_prepped(directories, molecule_name) == True:
            pass
         if self.is_poly_prepped(directories, molecule_name) == False:
            print("Please prepare prepin files for so polymers can be generated.")
            return(None)

         if number_of_units < 3:
             print("A minimum number of 3 units is required to construct the polymer.")
             return(None)
        
         molecule_dir = os.path.join(directories.molecules_dir, molecule_name)
         cd_command = "cd " + molecule_dir
         print(cd_command)

         molecule_dir = os.path.join(directories.molecules_dir, molecule_name)

         try:
             os.chdir(molecule_dir)
             print("Current directory:", os.getcwd())
         except Exception as e:
             print("Exception:", e)
         
         #try:
          #     result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
           #    if result.returncode == 0:
            #       # Command executed successfully
             #     print("Output:", result.stdout)
              # else:
               #    # Command failed, print error message
                #   print("Error:", result.stderr)
          #except Exception as e:
               # Exception occurred during subprocess execution
              #print("Exception:", e)
         
         file_subtype = "_3_3_array_" + str(number_of_units) + "_polymer"
         #head_prepi_filepath = os.path.join(directories.molecules_dir, molecule_name, ("head_" + molecule_name + ".prepi"))
         head_prepi_filepath = "head_" + molecule_name + ".prepi"
    
         #mainchain_prepi_filepath = os.path.join(directories.molecules_dir, molecule_name, ("mainchain_" + molecule_name + ".prepi"))
         mainchain_prepi_filepath = "mainchain_" + molecule_name + ".prepi"
    
         #tail_prepi_filepath = os.path.join(directories.molecules_dir, molecule_name, ("tail_" + molecule_name + ".prepi"))
         tail_prepi_filepath = "tail_" + molecule_name + ".prepi"
    
         output_dir = os.path.join(directories.systems_dir, (molecule_name.split("_")[0] + file_subtype))
         if not os.path.exists(output_dir):
             os.makedirs(output_dir)

         # Also need the pdb_file for pairwise distance (box_size distance), and the pdb_file of the monomer (for x,y distance)
         #pdb_filepath = os.path.join(directories.molecules_dir, molecule_name, (molecule_name + ".pdb"))
         pdb_filepath = molecule_name + ".pdb"
    
         monomer_name = molecule_name.split("_")[0] + "_monomer"
         monomer_pdb_filepath = os.path.join(directories.molecules_dir, monomer_name, (monomer_name + ".pdb"))
         current_directory = os.getcwd()
         print("Current working directory:", current_directory)
         Mol = MolFromPDBFile(pdb_filepath)
         monomer = MolFromPDBFile(monomer_pdb_filepath)
         max_dist = self.max_pairwise_distance(monomer)
         print("monomer max dist is: ", max_dist)
         translate_distance = float((int(max_dist)*2))
         polymer_length = max_dist * number_of_units
    
         #max_dist_box = max_dist*(number_of_units) # +1 to ensure box is bigger than the polymer
         #print("polymer max dist is: ", max_dist_box)
         box_dist = int(float(polymer_length*3))

         polymer_name = molecule_name.split("_")[0] + "_" + str(number_of_units) + "_polymer"
         molecule_name_1 = polymer_name + "_1"
         molecule_name_2 = polymer_name + "_2"
         molecule_name_3 = polymer_name + "_3"
         molecule_name_4 = polymer_name + "_4"
         molecule_name_5 = polymer_name + "_5"
         molecule_name_6 = polymer_name + "_6"
         molecule_name_7 = polymer_name + "_7"
         molecule_name_8 = polymer_name + "_8"
         molecule_name_9 = polymer_name + "_9"

         translate_line_1 = "{0.0 0.0 0.0}"
         translate_line_2 = "{" + str(translate_distance) + " 0.0 0.0}"
         #translate_line_2 = "{0.0 0.0 " + str(translate_distance) + "}"
         translate_line_3 = "{" + str(-translate_distance) + " 0.0 0.0}"
         #translate_line_3 = "{0.0 0.0 " + str(-translate_distance) + "}"

         translate_line_4 = "{" + str(translate_distance) + " " + str(translate_distance) + " 0.0}"
         #translate_line_4 = "{0.0 " + str(translate_distance) + " " + str(translate_distance) + "}"
         translate_line_5 = "{" + str(-translate_distance) + " " + str(translate_distance) + " 0.0}"
         #translate_line_5 = "{0.0 " + str(translate_distance) + " " + str(-translate_distance) + "}"
         translate_line_6 = "{0.0 " + str(translate_distance) + " 0.0}"

         translate_line_7 = "{" + str(translate_distance) + " " + str(-translate_distance) + " 0.0}"
         #translate_line_7 = "{0.0 " + str(-translate_distance) + " " + str(translate_distance) + "}"
         
         translate_line_8 = "{" + str(-translate_distance) + " " + str(-translate_distance) + " 0.0}"
         #translate_line_8 = "{0.0 " + str(-translate_distance) + " " + str(-translate_distance) + "}"
         translate_line_9 = "{0.0 " + str(-translate_distance) + " 0.0}"

         combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + " " + molecule_name_5 + " " + molecule_name_6 + " " + molecule_name_7 + " " + molecule_name_8 + " " + molecule_name_9 + "}"

         base_mol_name = molecule_name.split("_")[0]
         #intleap_path = os.path.join(output_dir, (base_mol_name + file_subtype + ".intleap"))
         intleap_path = base_mol_name + file_subtype + ".intleap"
    
         prmtop_filepath =  os.path.join(output_dir, base_mol_name + file_subtype + "_" + str(box_dist) + ".prmtop")
         rst_filepath = os.path.join(output_dir, base_mol_name + file_subtype + "_" + str(box_dist) + ".rst7")

         unsolved_prmtop_filepath =  os.path.join(output_dir, "unsolved_" + base_mol_name + file_subtype + ".prmtop")
         unsolved_rst_filepath = os.path.join(output_dir, "unsolved_" + base_mol_name + file_subtype + ".rst7")
        
         three_three_array_pdb_filepath = os.path.join(output_dir, base_mol_name + file_subtype + "_" + str(box_dist) + ".pdb")
         unsolved_three_three_array_pdb_filepath = os.path.join(output_dir, "unsolved_" + base_mol_name + file_subtype + ".pdb")
    
         head_rescode, mainchain_rescode, tail_rescode = directories.retrieve_polymeric_rescodes("3HB_trimer")

         print("head_code ", head_rescode)
         print("tail code ", tail_rescode)
         print("mainchain_code ", mainchain_rescode)
         polymer_code = " ".join([head_rescode] + [mainchain_rescode] * (number_of_units - 2) + [tail_rescode])
         polymer_command = "{" + polymer_code + "}"

         file_content = f"""source leaprc.gaff
         source leaprc.water.fb3
         source leaprc.protein.ff14SB

         loadamberprep {head_prepi_filepath}
         loadamberprep {mainchain_prepi_filepath}
         loadamberprep {tail_prepi_filepath}

         list
            
         {molecule_name_1} = sequence {polymer_command}
         {molecule_name_2} = sequence {polymer_command}
         {molecule_name_3} = sequence {polymer_command}
         {molecule_name_4} = sequence {polymer_command}
         {molecule_name_5} = sequence {polymer_command}
         {molecule_name_6} = sequence {polymer_command}
         {molecule_name_7} = sequence {polymer_command}
         {molecule_name_8} = sequence {polymer_command}
         {molecule_name_9} = sequence {polymer_command}

         check {molecule_name_1}
         
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
         saveamberparm system {unsolved_prmtop_filepath} {unsolved_rst_filepath}
         savepdb system {unsolved_three_three_array_pdb_filepath}
    
         solvatebox system TIP3PBOX {box_dist}

         saveamberparm system {prmtop_filepath} {rst_filepath}
         savepdb system {three_three_array_pdb_filepath}
         quit
         """
    
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

         cd_command = "cd " + str(directories.main_dir)
         result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
         return()

    
        