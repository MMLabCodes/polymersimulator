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
    """
    A class for building and simulating molecular systems.

    Provides methods for converting SMILES strings to PDB files, running simulations,
    generating residue codes, calculating distances, and aligning molecules.

    Attributes:
        packmol_path (str): Path to the Packmol executable.
        manager (object): Manager object to handle file paths and directories.
    """
    packmol_path = "/home/dan/packmol-20.14.4-docs1/packmol-20.14.4-docs1/packmol"
    
    def __init__(self, manager):
        """
        Initializes the BuildSystems class with the specified manager object.

        Args:
            manager (object): The manager object for interacting with directories and files.
        """      
        self.manager = manager
    
    def SmilesToPDB(self, smiles_string, output_file):
        """
        Converts a SMILES string to a PDB file.

        Args:
            smiles_string (str): The SMILES string of the molecule.
            output_file (str): The name of the output PDB file.

        Returns:
            None: Writes the 3D structure of the molecule to a PDB file.
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

    def SmilesToPDB_GenResCode(self, smiles, name):
        """
        Converts a SMILES string to a PDB file and generates a residue code if not found.

        Args:
            smiles (str): The SMILES string of the molecule.
            name (str): The name of the molecule.

        Returns:
            None: Writes the PDB file and updates the residue code CSV.
        """
        forbidden_codes = ["AAA", "BBB", "CCC", "UNL", "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "SEC", "TRP", "TYR", "VAL"]
        # Load existing residue codes from the CSV file
        residue_codes = self.load_residue_codes(self.manager.residue_code_csv)
        # Check if the name or smiles is already in the database
        existing_entry = self.find_existing_entry(residue_codes, name, smiles)
        if existing_entry:
            # Use existing residue code
            residue_code = existing_entry[2]  # Assuming code is the third column
        else:
            # Generate a unique 3-letter residue code excluding forbidden codes
            residue_code = self.generate_unique_residue_code(residue_codes, forbidden_codes)
            # Update the CSV file with the new entry
            self.update_residue_codes_csv(name, smiles, residue_code, self.manager.residue_code_csv)
        # Replace default "UNL" codes in the PDB content - NOTE, "UNL" is also a forbidden code as it is the default code.
        pdb_filepath = os.path.join(self.manager.pdb_file_dir, (name + ".pdb"))
        self.SmilesToPDB(smiles, pdb_filepath)    
        with open(pdb_filepath, "r") as pdb_file:
            lines = pdb_file.readlines()     
            for i, line in enumerate(lines):
                lines[i] = line.replace(" UNL", f" {residue_code}")   
        with open(pdb_filepath, "w") as pdb_file:              
            pdb_file.writelines(lines)
        return(pdb_filepath)

    def PDBToXYZ(self, molecule_name):
        # Get pdb filepath
        pdb_file = self.manager.load_pdb_filepath(molecule_name)

        # Make xyz filepath
        xyz_file = pdb_file.replace(".pdb", ".xyz")
    
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "xyz")

        mol = openbabel.OBMol()
        if not obConversion.ReadFile(mol, pdb_file):
            raise ValueError(f"Failed to read PDB file: {pdb_file}")

        if not obConversion.WriteFile(mol, xyz_file):
            raise ValueError(f"Failed to write XYZ file: {xyz_file}")

        print(f"Conversion successful: {pdb_file} → {xyz_file}")
        return(xyz_file)
    
    def load_residue_codes(self, residue_code_csv):
        """
        Loads existing residue codes from a CSV file.

        Args:
            residue_code_csv (str): The path to the residue code CSV file.

        Returns:
            list: List of residue codes from the CSV file.
        """
        # Load existing residue codes from the CSV file
        residue_codes = []	
        with open(residue_code_csv, "r") as csv_file:
            reader = csv.reader(csv_file)
            residue_codes = [row for row in reader]
        return residue_codes

    def generate_unique_residue_code(self, residue_codes, forbidden_codes=None):
        """
        Generates a unique 3-letter residue code not in the forbidden or existing codes.

        Args:
            residue_codes (list): List of existing residue codes.
            forbidden_codes (list, optional): List of forbidden residue codes.

        Returns:
            str: A unique 3-letter residue code.
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

    def GenRescode_4_PolyUnits(self, name):
        """
        Generates residue codes for polymeric units (e.g., trimers) and updates the CSV.

        Args:
            name (str): The name of the trimer molecule.

        Returns:
            None: Updates the residue code CSV.
        """
        if "trimer" not in name:
            print("Polymeric unit generation requires trimers. Please consult the build systems guide for information on how to do this")
            return()
    
        forbidden_codes = ["AAA", "BBB", "CCC", "UNL", "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "SEC", "TRP", "TYR", "VAL"]
   
        # Load existing residue codes from the CSV file
        residue_codes = self.load_residue_codes(self.manager.residue_code_csv)
    
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
    
        self.update_residue_codes_csv(head_name, head_smiles, head_code, self.manager.residue_code_csv)
        self.update_residue_codes_csv(mainchain_name, mainchain_smiles, mainchain_code, self.manager.residue_code_csv)
        self.update_residue_codes_csv(tail_name, tail_smiles, tail_code, self.manager.residue_code_csv)
        print("Head code assigned: ", head_code)
        print("Mainchain code assigned: ", mainchain_code)
        print("Tail code assigned: ", tail_code)
    
    def update_residue_codes_csv(self, name, smiles, residue_code, residue_code_csv):
        """
        Updates the CSV file with a new residue code entry.

        Args:
            name (str): The name of the molecule.
            smiles (str): The SMILES representation of the molecule.
            residue_code (str): The residue code to be added.
            residue_code_csv (str): Path to the CSV file for residue codes.

        Returns:
            None: Appends the new residue code to the CSV file.
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
        Finds an existing residue code entry by molecule name or SMILES string.

        Args:
            residue_codes (list): The list of existing residue codes.
            name (str): The name of the molecule.
            smiles (str, optional): The SMILES string of the molecule.

        Returns:
            list or None: The matching entry if found, otherwise None.
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


    def run_packmol(self, input_file_name):
        """
        Runs the Packmol simulation with the specified input file.

        Args:
            input_file_name (str): The name of the Packmol input file.

        Returns:
            None: Executes Packmol simulation using the provided input file.
        """
        if os.path.exists(self.packmol_path):
            print(f"Packmol executable exists at '{self.packmol_path}'.")
        else:
            print(f"Packmol executable not found at '{self.packmol_path}'.")
            print("Please update the class variable for the packmol path using 'update_packmol_path(new_packmol_path)'")
            print("")
            print("Call packmol_help() for more information")
            return()
        packmol_filepath = self.manager.load_pckml_filepath(input_file_name)
        if os.path.exists(packmol_filepath):
            print(f"Packmol input file exists at '{packmol_filepath}'.")
        else:
            print("Packmol input file not found.")
            print("")
            return()
        system_dir = os.path.join(self.manager.systems_dir, input_file_name)
        packmol_command = self.packmol_path + " < " + packmol_filepath
        result = subprocess.run(packmol_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return()

    @classmethod
    def update_packmol_path(cls, new_packmol_path):
        """
        Updates the path to the Packmol executable.

        Args:
            new_packmol_path (str): The new path to the Packmol executable.

        Returns:
            None: Updates the class-level variable for the Packmol path.
        """
        if os.path.exists(new_packmol_path):
            self.packmol_path = new_packmol_path
            print(f"Packmol path updated to '{self.packmol_path}'.")
        else:
            print(f"The provided path '{new_packmol_path}' does not exist. Please provide a valid path.")

    @staticmethod
    def get_xyz_dists(input_file=None):
        """
        Calculates the maximum distance between the largest and smallest XYZ coordinates in a file.

        Args:
            input_file (str): Path to a '.pdb' or '.xyz' file.

        Returns:
            tuple: Maximum distances along the x, y, and z axes.
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
        """
        Displays help information for the `get_xyz_dists` method.

        Returns:
            None: Prints out the docstring of `get_xyz_dists` method.
        """
        print(cls.get_xyz_dists.__doc__)
        print("This method calculates the maximum distance between the largest and smallest xyz coordinates")
        print("from a PDB or XYZ file.")
           
    @staticmethod
    def align_molecule(input_pdb, axis_to_align=None):
        """
        Aligns the molecule structure to a specified axis (X, Y, or Z).

        Args:
            input_file (str): Path to a PDB file containing the molecule.
            axis (str): The axis to align the molecule to ('X', 'Y', or 'Z').

        Returns:
            None: The aligned structure is saved with a new file name.
        """
        import MDAnalysis as mda
        import numpy as np
        from scipy.spatial.transform import Rotation as R
            
        u = mda.Universe(input_pdb)
            
        # Select the atom group for alignment (you can specify your selection if needed)
        ag = u.atoms

        # Compute the principal axes using the mass-weighted inertia tensor
        principal_axes = ag.principal_axes()

        # Define the target axis, which is along the z-axis ([0, 0, 1])
        z_axis = np.array([0, 0, 1])
        x_axis = np.array([1, 0, 0])
        y_axis = np.array([0, 1, 0])
        
        # Get the principal axis that we want to align with the z-axis (typically the first one, corresponding to the largest eigenvalue)
        principal_axis_to_align = principal_axes[0]

        # Compute the rotation needed to align the principal axis with the z-axis
        if axis_to_align == "Z":
            axis = z_axis
        elif axis_to_align == "X":
            axis = x_axis
        elif axis_to_align == "Y":
            axis = y_axis
        else:
            print("Please provide arugment as: 'X' OR 'Y' OR 'Z'")
            
        rotation, _ = R.align_vectors([principal_axis_to_align], [axis])

        # Apply the rotation to the entire molecule
        ag.positions = rotation.apply(ag.positions)

        # Save the aligned molecule to a new PDB file
        u.atoms.write(input_pdb)

        return(input_pdb)
            
            

class BuildAmberSystems(BuildSystems):
    
    error_param = "Molecule not parametrized. Please parametrize pdb_file."
    
    def __init__(self, manager):
        self.manager = manager
    
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
    
    def parameterize_mol(self, molecule_name=None, forcefield=None):
        if  molecule_name == None:
            print("Please provide 1 argument as follows: parametrize_mol(molecule_name)")
            print("Directories: A python object generated with the PolymerSimulatorDirs(filepath) method imported from sw_directories")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            return(None)

        if forcefield == None:
            forcefield = "GAFF"
        else:
            forcefield = forcefield
        # Create a new directory for param files for the molecule and copy pdb there
        pdb_filepath = os.path.join(self.manager.pdb_file_dir, (molecule_name + ".pdb"))
        self.mod_pdb_file(pdb_filepath)
        param_mol_dir = os.path.join(self.manager.molecules_dir, molecule_name)
        if not os.path.exists(param_mol_dir):
            os.makedirs(param_mol_dir, exist_ok=True)
        shutil.copy2(pdb_filepath, param_mol_dir)
        os.remove(pdb_filepath)
        
        # Specify paths for tleap
        pdb_filepath = os.path.join(param_mol_dir, (molecule_name + ".pdb"))
        mol2_filepath = os.path.join(param_mol_dir, (molecule_name + ".mol2"))
        antechamber_command = "antechamber -i " + pdb_filepath + " -fi pdb -o " + mol2_filepath + " -fo mol2 -c bcc -s 2"       
        subprocess.run(antechamber_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
        frcmod_filepath = os.path.join(param_mol_dir, (molecule_name + ".frcmod"))
        parmchk_command = "parmchk2 -i " + mol2_filepath + " -f mol2 -o " + frcmod_filepath      
        subprocess.run(parmchk_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        lib_filepath = os.path.join(param_mol_dir, (molecule_name + ".lib"))
        intleap_filepath = os.path.join(param_mol_dir, (molecule_name + ".intleap"))
        
        file_content = f"""source leaprc.protein.ff14SB
        source leaprc.{forcefield}
        source leaprc.water.fb3
        {molecule_name} = loadmol2 {mol2_filepath}
        saveoff {molecule_name} {lib_filepath}
        quit
        """
        with open(intleap_filepath, 'w') as file:
            file.write(file_content)
        leap_command = "tleap -f " + intleap_filepath
        subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.gen_prepin_file_sing_mol(molecule_name)

    def gen_ac_file(self, molecule_name=None):
        # This function isn't required so much - however it will generate .ac files 
        if molecule_name == None:
            print("Please provide 1 argument as follows: gen_ac_file(molecule_name)")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            return(None)
        mol2_name = os.path.join(self.manager.molecules_dir, molecule_name, molecule_name) + ".mol2"
        ac_name = os.path.join(self.manager.molecules_dir, molecule_name, molecule_name) + ".ac"
        antechamber_command = "antechamber -fi mol2 -fo ac -i " + mol2_name + " -o " + ac_name + " -c bcc -s 2"
        subprocess.run(antechamber_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return(ac_name)

    def gen_prepin_files(self, molecule_name):
        # NOTE: this function is for generating prepin files of polymers, for single molecules, please use 'gen_prepin_files_single_mol'
        if molecule_name == None:
            print("Please provide 1 argument as follows: gen_prepin_files(molecule_name)")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            return(None)
        
        # Find the path of the molecule directory (where all the required files are required)
        molecule_dir = os.path.join(self.manager.molecules_dir, molecule_name)
        if os.path.exists(molecule_dir):
            pass
        else:
            print("Please parameterize molecule.")
            return()  
    
        # Define file names required for prepgen
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
        head_rescode, mainchain_rescode, tail_rescode = self.manager.retrieve_polymeric_rescodes(molecule_name)
        print(head_rescode, mainchain_rescode, tail_rescode)
        if head_rescode == None or mainchain_rescode == None or tail_rescode == None:
            print("Residue codes for polymeric units not generated. Please generate them with 'build.GenRescode_4_PolyUnits('ethane')'")
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
        os.chdir(self.manager.main_dir)
        return(None)

    def gen_prepin_file_sing_mol(self, molecule_name):
        if molecule_name == None:
            print("Please provide 1 argument as follows: gen_prepin_files(molecule_name)")
            print("Molecule name: A string of the molecule name, i.e. 'Ethane'")
            return(None)
            
        # Find the path of the molecule directory (where all the required files are required)
        molecule_dir = os.path.join(self.manager.molecules_dir, molecule_name)
        if os.path.exists(molecule_dir):
            pass
        else:
            print("Please parameterize molecule.")
            return()  
    
        # Define file names required for prepgen
        prepi_file = os.path.join(molecule_dir, molecule_name + ".prepi")
        pdb_file = os.path.join(molecule_dir, molecule_name + ".pdb")
        frcmod_file = os.path.join(molecule_dir, molecule_name + ".frcmod")

        antechamber_command = f"antechamber -i {pdb_file} -fi pdb -o {prepi_file} -fo prepi -c bcc -s 2"
        try:
            result = subprocess.run(antechamber_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
                # Command executed successfully
                print("Output:", result.stdout)
            else:
                # Command failed, print error message
                print("Error:", result.stderr)
        except Exception as e:
            # Exception occurred during subprocess execution
            print("Exception:", e)

        parmchk_command = f"parmchk2 -i {prepi_file} -f prepi -o {frcmod_file}"
        try:
            result = subprocess.run(antechamber_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
                # Command executed successfully
                print("Output:", result.stdout)
            else:
                # Command failed, print error message
                print("Error:", result.stderr)
        except Exception as e:
            # Exception occurred during subprocess execution
            print("Exception:", e)
        return(None)
    
    def is_mol_parametrized(self,  molecule_name):
        param_mol_dir = os.path.join(self.manager.molecules_dir, molecule_name)
        molecule_dir = os.path.join(self.manager.pdb_file_dir, molecule_name)
        pdb_filepath = os.path.join(param_mol_dir, (molecule_name + ".pdb"))
        mol2_filepath = os.path.join(param_mol_dir, (molecule_name + ".mol2"))
        files_exist = os.path.exists(pdb_filepath) and os.path.exists(mol2_filepath)
        return(files_exist) # This will be true or false

    def is_poly_prepped(self, molecule_name):
        poly_prepped_dir = os.path.join(self.manager.molecules_dir, molecule_name)
        head_prepi_filepath = os.path.join(poly_prepped_dir, ("head_" + molecule_name + ".prepi"))
        mainchain_prepi_filepath = os.path.join(poly_prepped_dir, ("mainchain_" + molecule_name + ".prepi"))
        tail_prepi_filepath = os.path.join(poly_prepped_dir, ("tail_" + molecule_name + ".prepi"))
        files_exist = os.path.exists(head_prepi_filepath) and os.path.exists(mainchain_prepi_filepath) and os.path.exists(tail_prepi_filepath)
        return(files_exist) # This will be true or false

    def pdb_2_mol2(self, pdb_file):
        mol2_file = (pdb_file.split(".")[0]) + ".mol2"
        babel_command = "obabel -ipdb " + pdb_file + " -omol2 -O " + mol2_file + " --partialcharge gasteiger"
        try:
            result = subprocess.run(babel_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
               # Command executed successfully
               pass
            else:
                # Command failed, print error message
                print("Error:", result.stderr)
        except Exception as e:
            # Exception occurred during subprocess execution
            print("Exception:", e)  
        return(None)

    def relax_molecule(self, mol, pdb_file, forcefield=None, max_iters=200):
        """
        Relax the molecule using the specified force field.

        This isn't useful for polymers as wiggly conformations can be generated, this will also fail when there are very large molecules with many degrees of freedom
        It is advised to used this function for small molecules.
        """
        AllChem.EmbedMolecule(mol)
    
        # Optimize the geometry
        if forcefield == 'UFF':
            AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
        elif forcefield == 'MMFF':
            AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
        else:
            raise ValueError(f"Unknown force field: {forcefield}")
    
        # Write the relaxed molecule to the output PDB file
        with Chem.PDBWriter(pdb_file) as writer:
            writer.write(mol)
        
        return mol

    def extract_coordinates_from_pdb(self, pdb_file):
        """
        Extracts coordinates from a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.

        Returns:
            list of tuples: A list of (x, y, z) coordinates for each atom.
        """
        coordinates = []
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates.append((x, y, z))
        return coordinates

    def replace_coordinates_in_pdb(self, original_pdb_file, pdb_file_with_new_coords):
        """
        Replaces the coordinates in the original PDB file with new coordinates and writes to a new file.

        Args:
            original_pdb_file (str): Path to the original PDB file.
            pdb_file_with_new_coords (str): Path to the PDB file with the new coordinates.
        """
        new_coordinates = self.extract_coordinates_from_pdb(pdb_file_with_new_coords)
    
        with open(original_pdb_file, 'r') as infile:
            lines = infile.readlines()

        # Modify the lines with new coordinates
        updated_lines = []
        coord_index = 0
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if coord_index < len(new_coordinates):
                    new_x = f"{new_coordinates[coord_index][0]:8.3f}"
                    new_y = f"{new_coordinates[coord_index][1]:8.3f}"
                    new_z = f"{new_coordinates[coord_index][2]:8.3f}"
                    new_line = f"{line[:30]}{new_x:>8}{new_y:>8}{new_z:>8}{line[54:]}"
                    coord_index += 1
                    #print(new_line)
                else:
                    new_line = line
                updated_lines.append(new_line)
            else:
                updated_lines.append(line)
       
        # Write the updated content back to the original file
        with open(original_pdb_file, 'w') as outfile:
            outfile.writelines(updated_lines)
    
    def gen_polymer_pdb_and_params(self, molecule_name, number_of_units, box_radius=None, infinite=None):
        """
        Generates a polymer PDB file using `tleap` based on the specified molecule and number of units.

        Args:
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
        if box_radius == None:
            box_radius = 10.0
        else:
            box_radius = box_radius
            if type(box_radius) is float:
                pass
            else:
                print("Plase pass the box radius as a float")
                print("Example: gen_amber_params_sing_poly('base_trimer_name', 'polymer_name', 20.0)")
                return()
                
        molecule_dir = os.path.join(self.manager.molecules_dir, molecule_name)
        try:
            os.chdir(molecule_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
        
        file_subtype = "_" + str(number_of_units) + "_polymer"
        head_prepi_filepath = "head_" + molecule_name + ".prepi"
        tail_prepi_filepath = "tail_" + molecule_name + ".prepi"
        mainchain_prepi_filepath = "mainchain_" + molecule_name + ".prepi"

        output_dir = os.path.join(self.manager.systems_dir, (molecule_name.split("_")[0] + file_subtype))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        polymer_name = molecule_name.split("_")[0] + "_" + str(number_of_units) + "_polymer"
        if infinite is not None:
            polymer_name = polymer_name + "_infinite"

        intleap_path = polymer_name + ".intleap"
        prmtop_filepath = os.path.join(output_dir, polymer_name + ".prmtop")
        rst_filepath = os.path.join(output_dir, polymer_name + ".rst7")
        pdb_filepath = os.path.join(output_dir, polymer_name + ".pdb")

        head_rescode, mainchain_rescode, tail_rescode = self.manager.retrieve_polymeric_rescodes(molecule_name)

        if infinite is not None:
            polymer_code = " ".join([mainchain_rescode]*number_of_units) 
        else:
            polymer_code = " ".join([head_rescode] + [mainchain_rescode] * (number_of_units - 2) + [tail_rescode])
        polymer_command = "{" + polymer_code + "}"

        file_content = f"""loadamberprep {head_prepi_filepath}
             loadamberprep {mainchain_prepi_filepath}
             loadamberprep {tail_prepi_filepath}

             list

             polymer = sequence {polymer_command}
             setBox polymer centers {box_radius}
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
            os.chdir(self.manager.main_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
        return(pdb_filepath, prmtop_filepath, rst_filepath)
        
    def gen_copolymer_pdb_and_params(self, pattern, base_trimers, number_of_units, box_radius=None):
        # Check base_trimers were passed as a list
        if isinstance(base_trimers, list):
            pass
        else:
            print(f"Arguments should be passed the function as follows.")
            print("Pattern (string): 'AB'")
            print("Base trimers (list): ['3HB_trimer', '3HV_trimer']")
            print("Number of units (int): 10")
            return

        if box_radius == None:
            box_radius = 0.0
        else:
            box_radius = box_radius
            if type(box_radius) is float:
                pass
            else:
                print("Plase pass the box radius as a float")
                print("Example: gen_copolymer_pdb_and_params('AB', ['base_trimer_1', 'base_trimer_2'], 10, 20.0)")
                return()
            

        # Create the output dircotry and name the system
        file_subtype = "_copolymer_" + pattern
        copolymer_name = "_".join([item.split("_")[0] for item in base_trimers]) + "_" + str(number_of_units) + file_subtype 
        print(copolymer_name)  
    
        # If only a copolymer of 2 polymers (and not 3) append the first polymers to the list of base trimers
        if len(base_trimers) == 3:
            pass
        if len(base_trimers) == 2:
            # Need to append something else to the list so the function works, but this will not be needed (it loads 3 polymers every time and the code will break if not, easiest thing is to load the third to be the same as one of the others)
            base_trimers.append(base_trimers[0])

        # Get paths to prepi files for head, mainchain and tail of the desired trimers
        base_trimer_dirs = [os.path.join(self.manager.molecules_dir, base_trimers[0]), os.path.join(self.manager.molecules_dir, base_trimers[1]), os.path.join(self.manager.molecules_dir, base_trimers[2])]
        head_prepi_filenames = [f"head_{base_trimers[0]}.prepi", f"head_{base_trimers[1]}.prepi", f"head_{base_trimers[2]}.prepi"]
        main_prepi_filenames = [f"mainchain_{base_trimers[0]}.prepi", f"mainchain_{base_trimers[1]}.prepi", f"mainchain_{base_trimers[2]}.prepi"]
        tail_prepi_filenames = [f"tail_{base_trimers[0]}.prepi", f"tail_{base_trimers[1]}.prepi", f"tail_{base_trimers[2]}.prepi"]
        head_prepi_filepaths = [os.path.join(base_trimer_dirs[0], head_prepi_filenames[0]), os.path.join(base_trimer_dirs[1], head_prepi_filenames[1]), os.path.join(base_trimer_dirs[2], head_prepi_filenames[2])]
        main_prepi_filepaths = [os.path.join(base_trimer_dirs[0], main_prepi_filenames[0]), os.path.join(base_trimer_dirs[1], main_prepi_filenames[1]), os.path.join(base_trimer_dirs[2], main_prepi_filenames[2])]
        tail_prepi_filepaths = [os.path.join(base_trimer_dirs[0], tail_prepi_filenames[0]), os.path.join(base_trimer_dirs[1], tail_prepi_filenames[1]), os.path.join(base_trimer_dirs[2], tail_prepi_filenames[2])]

        # Get paths to the .frcmod files
        base_trimer_frcmods = [os.path.join(base_trimer_dirs[0], f"{base_trimers[0]}.frcmod"), os.path.join(base_trimer_dirs[1], f"{base_trimers[1]}.frcmod"), os.path.join(base_trimer_dirs[2], f"{base_trimers[2]}.frcmod")]
 
        output_dir = os.path.join(self.manager.systems_dir, copolymer_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Create files for use with tleap
        intleap_path = os.path.join(output_dir, copolymer_name + ".intleap")
        rst_path = os.path.join(output_dir, copolymer_name + ".rst7")
        prmtop_path = os.path.join(output_dir, copolymer_name + ".prmtop")
        pdb_path = os.path.join(output_dir, copolymer_name + ".pdb")

        # Retrieve residue codes for each polymeric unit
        A_head_rescode, A_mainchain_rescode, A_tail_rescode = self.manager.retrieve_polymeric_rescodes(base_trimers[0])
        B_head_rescode, B_mainchain_rescode, B_tail_rescode = self.manager.retrieve_polymeric_rescodes(base_trimers[1])
        C_head_rescode, C_mainchain_rescode, C_tail_rescode = self.manager.retrieve_polymeric_rescodes(base_trimers[2])
    
        # Store residue codes in a dictionary
        rescode_dict = {
            "A": {
                "head": A_head_rescode,
                "mainchain": A_mainchain_rescode,
                "tail": A_tail_rescode
            },
            "B": {
                "head": B_head_rescode,
                "mainchain": B_mainchain_rescode,
                "tail": B_tail_rescode
            },
            "C": {
                "head": C_head_rescode,
                "mainchain": C_mainchain_rescode,
                "tail": C_tail_rescode
            }
        }

        # This section of code makes a pattern for the copolymer (i.e. if a decamer in 'AB' configuration, the pattern is 'ABABABABAB')
        if len(pattern) == number_of_units:
            making_pattern = pattern
        elif (number_of_units / len(pattern)).is_integer():
            making_pattern = pattern * int(number_of_units / len(pattern))
        else:
            print("Pattern and number of units incompatible. Please revise your pattern and retry this function")

        # Parse the pattern and attach the specific residue codes 
        copolymer_pattern = ""
        for i in range(number_of_units):
            if i == 0:
                copolymer_pattern = rescode_dict[making_pattern[i]]["head"]
            elif i != number_of_units-1:
                copolymer_pattern = copolymer_pattern + " " + rescode_dict[making_pattern[i]]["mainchain"]
            else:
                copolymer_pattern = copolymer_pattern + " " + rescode_dict[making_pattern[i]]["tail"]

        copolymer_command = "{" + copolymer_pattern + "}"
        print(copolymer_command)

        file_content = f"""source leaprc.gaff
            source leaprc.water.fb3
            source leaprc.protein.ff14SB

            loadamberprep {head_prepi_filepaths[0]}
            loadamberprep {head_prepi_filepaths[1]}
            loadamberprep {head_prepi_filepaths[2]}
            loadamberprep {main_prepi_filepaths[0]}
            loadamberprep {main_prepi_filepaths[1]}
            loadamberprep {main_prepi_filepaths[2]}
            loadamberprep {tail_prepi_filepaths[0]}
            loadamberprep {tail_prepi_filepaths[1]}
            loadamberprep {tail_prepi_filepaths[2]}

            loadamberparams {base_trimer_frcmods[0]}
            loadamberparams {base_trimer_frcmods[1]}
            loadamberparams {base_trimer_frcmods[2]}

            list

            polymer = sequence {copolymer_command}
            setBox polymer centers {box_radius}
            saveamberparm polymer {prmtop_path} {rst_path}
            savepdb polymer {pdb_path}
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

        return(copolymer_name)

    def gen_amber_params_sing_mol_solvated(self, molecule_name, buffer=None):
        if buffer == None:
            buffer = "10"
        else:
            buffer = str(buffer)
        
        molecule_dir = os.path.join(self.manager.molecules_dir, molecule_name)

        try:
            os.chdir(molecule_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
            return(None)

        molecule_prepi_filepath = molecule_name + ".prepi"
        molecule_frcmod_filepath = molecule_name + ".frcmod"
    
        output_dir = os.path.join(self.manager.systems_dir, (molecule_name + "_wat_solv"))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        molecule_pdb_filepath = os.path.join(self.manager.molecules_dir, molecule_name, (molecule_name + ".pdb"))
        x, y, z = self.get_xyz_dists(molecule_pdb_filepath)

        file_subtype = "_wat_solv"
        filename = molecule_name + file_subtype + f"_{buffer}"
        intleap_path = filename + ".intleap"

        prmtop_filepath = os.path.join(output_dir, filename + ".prmtop")
        rst_filepath = os.path.join(output_dir, filename + ".rst7")
        pdb_output = os.path.join(output_dir, filename + ".pdb")

        frcmod_filepath = os.path.join(self.manager.molecules_dir, molecule_name, molecule_name + ".frcmod")
        lib_filepath = os.path.join(self.manager.molecules_dir, molecule_name, molecule_name + ".lib")
        box_dist_line = "solvatebox " + molecule_name + " TIP3PBOX {" + buffer + "}"
               
        file_content = f"""source leaprc.gaff
            source leaprc.water.fb3
            loadamberparams {frcmod_filepath}
            loadoff {lib_filepath}

            solvatebox {molecule_name} TIP3PBOX {buffer}

            saveamberparm {molecule_name} {prmtop_filepath} {rst_filepath}
            savepdb {molecule_name} {pdb_output}
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

        os.chdir(self.manager.main_dir)
        return(filename)

    def gen_amber_params_poly_solvated(self, base_molecule_name, polymer_name, buffer=None, water_model=None):

        pdb_file = self.manager.load_pdb_filepath(polymer_name)
        base_pdb_file = self.manager.load_pdb_filepath(base_molecule_name)
         
        base_pdb_dir = os.path.join(self.manager.molecules_dir, base_molecule_name)

        head_prepi_filepath = os.path.join(base_pdb_dir, "head_" + base_molecule_name + ".prepi")
        mainchain_prepi_filepath = os.path.join(base_pdb_dir, "mainchain_" + base_molecule_name + ".prepi")
        tail_prepi_filepath = os.path.join(base_pdb_dir, "tail_" + base_molecule_name + ".prepi")
        frcmod_filepath = os.path.join(base_pdb_dir, base_molecule_name + ".frcmod")
       
        if buffer == None:
            buffer = "10"
        else:
            buffer = str(buffer)

        if water_model == None:
            water_model = "TIP3PBOX"
        else:
            water_model = water_model
        
        output_dir = os.path.join(self.manager.systems_dir, (polymer_name + "_wat_solv"))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        file_subtype = "_wat_solv"
        filename = polymer_name + file_subtype + f"_{buffer}"
        intleap_path = filename + ".intleap"

        prmtop_filepath = os.path.join(output_dir, filename + ".prmtop")
        rst_filepath = os.path.join(output_dir, filename + ".rst7")
        pdb_output = os.path.join(output_dir, filename + ".pdb")
        
        file_content = f"""source leaprc.gaff
            source leaprc.water.fb3
            loadamberprep {head_prepi_filepath}
            loadamberprep {mainchain_prepi_filepath}
            loadamberprep {tail_prepi_filepath}
            loadamberparams {frcmod_filepath}

            polymer = loadpdb {pdb_file}
            
            solvatebox polymer {water_model} {buffer}

            saveamberparm polymer {prmtop_filepath} {rst_filepath}
            savepdb polymer {pdb_output}
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
        return(filename)
    
    def generate_3_3_polymer_array_crystal(self, base_molecule_name=None, molecule_name=None):
        # Thsis function builds arrays of polymers using the pre generated pdb files
        if molecule_name == None or base_molecule_name == None:
         # UPDATE THIS
            print("Please provide 2 arguments as follows: build_3_3_polymer_array(base_molecule_nmae, molecule_name)")
            print("Base polymer name: A string of the polymer name, i.e. '3HB_trimer'")
            print("Polymer name: A string of the polymer name, i.e. '3HB_10_polymer'")
            
            return(None)        

        pdb_file = self.manager.load_pdb_filepath(molecule_name)
        base_pdb_file = self.manager.load_pdb_filepath(base_molecule_name)
         
        molecule_dir = os.path.join(self.manager.molecules_dir, base_molecule_name)
        cd_command = "cd " + molecule_dir

        try:
            os.chdir(molecule_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
             
        head_prepi_filepath = "head_" + base_molecule_name + ".prepi"
        mainchain_prepi_filepath = "mainchain_" + base_molecule_name + ".prepi"
        tail_prepi_filepath = "tail_" + base_molecule_name + ".prepi"

        file_subtype = "_3_3_array_crystal"
        system_name = molecule_name + file_subtype
        output_dir = os.path.join(self.manager.systems_dir, system_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        x, y, z = self.get_xyz_dists(base_pdb_file)

        translate_distance = int((max(x, y))/2) # Removed z as they should not overlap in this distance
        z_trans = (z/2)+2
        y_trans = (y/2)+2

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
        translate_line_2 = "{0.0 0.0 " + str(z_trans) + " }"
        translate_line_3 = "{0.0 0.0 " + str(-z_trans) + " }"

        translate_line_4 = "{0.0 " + str(y_trans) + " " + str(z_trans) + "}"
        translate_line_5 = "{0.0 " + str(y_trans) + " 0.0}"
        translate_line_6 = "{0.0 " + str(y_trans) + " " + str(-z_trans) + "}"

        translate_line_7 = "{0.0 " + str(-y_trans) + " " + str(z_trans) + "}"
        translate_line_8 = "{0.0 " + str(-y_trans) + " 0.0}"
        translate_line_9 = "{0.0 " + str(-y_trans) + " " + str(-z_trans) + "}"

        combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + " " + molecule_name_5 + " " + molecule_name_6 + " " + molecule_name_7 + " " + molecule_name_8 + " " + molecule_name_9 + "}"

        base_mol_name = molecule_name.split("_")[0]
        intleap_path = system_name + ".intleap"

        prmtop_filepath =  os.path.join(output_dir, system_name + ".prmtop")
        rst_filepath = os.path.join(output_dir, system_name + ".rst7")
     
        three_three_array_pdb_filepath = os.path.join(output_dir, system_name + ".pdb")
 
        head_rescode, mainchain_rescode, tail_rescode = self.manager.retrieve_polymeric_rescodes(base_molecule_name)

        file_content = f"""source leaprc.gaff
        source leaprc.water.fb3
        source leaprc.protein.ff14SB

        loadamberprep {head_prepi_filepath}
        loadamberprep {mainchain_prepi_filepath}
        loadamberprep {tail_prepi_filepath}

        list
            
        {molecule_name_1} = loadpdb {pdb_file}
        {molecule_name_2} = loadpdb {pdb_file}
        {molecule_name_3} = loadpdb {pdb_file}
        {molecule_name_4} = loadpdb {pdb_file}
        {molecule_name_5} = loadpdb {pdb_file}
        {molecule_name_6} = loadpdb {pdb_file}
        {molecule_name_7} = loadpdb {pdb_file}
        {molecule_name_8} = loadpdb {pdb_file}
        {molecule_name_9} = loadpdb {pdb_file}

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
        setBox system vdw 0.0
        saveamberparm system {prmtop_filepath} {rst_filepath}
        savepdb system {three_three_array_pdb_filepath}
    
        quit
        """
    
        with open(intleap_path, 'w') as file:
            file.write(file_content)
            
        leap_command = "tleap -f " + intleap_path
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

        cd_command = "cd " + str(self.manager.main_dir)
        result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return(system_name)
    
    def generate_3_3_polymer_array_random(self, base_molecule_name, molecule_name):
        """
        Generates a 3x3 array of polymers using Packmol and adds 'TER' records after each polymer chain.

        This function creates a 3x3 array of polymer molecules based on a given input molecule. It uses Packmol 
        to arrange the molecules in the specified array and adds 'TER' records after each polymer chain to ensure 
        proper formatting for subsequent simulations.

        Parameters:
        dirs (object): An object that provides methods for directory and file management.
        molecule_name (str): The name of the molecule to be used for creating the polymer array. i.e. "3HB_10_polymer"
        base_molecule_name (str): The base name of the molecule used to identify the polymeric residue codes. i.e. "3HB_trimer"

        Returns:
        None

        The function performs the following steps:
        1. Defines the file subtype and constructs the PDB file path.
        2. Creates an output directory if it doesn't exist.
        3. Constructs the output PDB file name for the unsolved array.
        4. Calculates the dimensions of the box required to contain the 3x3 array of polymers.
        5. Generates the Packmol input file content with the specified parameters.
        6. Writes the Packmol input file to the output directory.
        7. Executes the Packmol command to generate the unsolved array PDB file.
        8. Adds 'TER' records to the generated PDB file using the `add_ter_to_pckml_result` function.

        Example usage:
            generate_3_3_polymer_array_pckml(dirs, 'polymer_molecule', 'base_polymer')
        """
        file_subtype = "_3_3_array_random"
        pdb_file = self.manager.load_pdb_filepath(molecule_name)
        system_name = molecule_name + file_subtype
        output_dir = os.path.join(self.manager.systems_dir, system_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        array_pdb_name = system_name + ".pdb"
        array_pdb = os.path.join(output_dir, array_pdb_name)
        x,y,z = self.get_xyz_dists(pdb_file)
        box_edge_len = int((max(x,y,z))*2)
        box_sizes = box_edge_len, box_edge_len, box_edge_len
        file_content=f"""tolerance 2.0
            output {array_pdb}
            filetype pdb
            structure {pdb_file}
            	number 9
                inside box 0. 0. 0. {box_sizes[0]}. {box_sizes[1]}. {box_sizes[2]}.
            end structure
            """
        packmol_input_file = system_name + ".pckml"
        packmol_input_filepath = os.path.join(output_dir, packmol_input_file)
        with open(packmol_input_filepath, 'w') as file:
            file.write(file_content)
        packmol_command = str(self.packmol_path) + " < " + packmol_input_filepath
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
        
        self.add_ter_to_pckml_result(array_pdb, base_molecule_name)
        return(system_name)

    def gen_amber_params_4_pckml_array(self, system_name, base_molecule_name):
        input_pdb = self.manager.load_pdb_filepath(system_name)
        
        head_prepi_filepath = os.path.join(self.manager.molecules_dir, base_molecule_name, ("head_" + base_molecule_name + ".prepi"))
        mainchain_prepi_filepath = os.path.join(self.manager.molecules_dir, base_molecule_name, ("mainchain_" + base_molecule_name + ".prepi"))
        tail_prepi_filepath = os.path.join(self.manager.molecules_dir, base_molecule_name, ("tail_" + base_molecule_name + ".prepi"))

        file_subtype = "_3_3_array_random"
        
        output_dir = os.path.join(self.manager.systems_dir, system_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        intleap_filepath = os.path.join(output_dir, system_name + ".intleap")
        prmtop_filepath =  os.path.join(output_dir, system_name + ".prmtop")
        rst_filepath = os.path.join(output_dir, system_name + ".rst7")
        head_rescode, mainchain_rescode, tail_rescode = self.manager.retrieve_polymeric_rescodes(base_molecule_name)

        file_content = f"""source leaprc.gaff

        loadamberprep {head_prepi_filepath}
        loadamberprep {mainchain_prepi_filepath}
        loadamberprep {tail_prepi_filepath}

        list
            
        system = loadpdb {input_pdb}

        check

        setBox system vdw 0.0
        saveamberparm system {prmtop_filepath} {rst_filepath}

        quit
        """
        with open(intleap_filepath, 'w') as file:
            file.write(file_content)
            
        leap_command = "tleap -f " + intleap_filepath
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

        return(system_name)
        
    def generate_polymer_3_3_array(self, base_molecule_name, molecule_name, method):
        if method == "crystal":
            syste_name = self.generate_3_3_polymer_array_crystal(base_molecule_name, molecule_name)
        if method == "random":
            system_name = self.generate_3_3_polymer_array_random(base_molecule_name, molecule_name)
            self.gen_amber_params_4_pckml_array(system_name, base_molecule_name)
        return(system_name)
            
    def generate_5_5_polymer_array_crystal(self, base_molecule_name=None, molecule_name=None, crystal_trans=None):
        # Thsis function builds arrays of polymers using the pre generated pdb files
        if molecule_name == None or base_molecule_name == None:
         # UPDATE THIS
            print("Please provide 2 arguments as follows: build_5_5_polymer_array(base_molecule_nmae, molecule_name)")
            print("Base polymer name: A string of the polymer name, i.e. '3HB_trimer'")
            print("Polymer name: A string of the polymer name, i.e. '3HB_10_polymer'")
            
            return(None)        

        pdb_file = self.manager.load_pdb_filepath(molecule_name)
        base_pdb_file = self.manager.load_pdb_filepath(base_molecule_name)
         
        molecule_dir = os.path.join(self.manager.molecules_dir, base_molecule_name)

        try:
            os.chdir(molecule_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
             
        head_prepi_filepath = "head_" + base_molecule_name + ".prepi"
        mainchain_prepi_filepath = "mainchain_" + base_molecule_name + ".prepi"
        tail_prepi_filepath = "tail_" + base_molecule_name + ".prepi"

        file_subtype = "_5_5_array_crystal"
        system_name = molecule_name + file_subtype
        output_dir = os.path.join(self.manager.systems_dir, system_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        x, y, z = self.get_xyz_dists(base_pdb_file)

        translate_distance = int((max(x, y))/2) # Removed z as they should not overlap in this distance
        if crystal_trans is None:
            z_trans = int((z/2)+2)
            y_trans = int((y/2)+2)
        else:
            z_trans = int((z/2)+crystal_trans)
            y_trans = int((y/2)+crystal_trans)

        molecule_name_1 = molecule_name + "_1"
        molecule_name_2 = molecule_name + "_2"
        molecule_name_3 = molecule_name + "_3"
        molecule_name_4 = molecule_name + "_4"
        molecule_name_5 = molecule_name + "_5"
        molecule_name_6 = molecule_name + "_6"
        molecule_name_7 = molecule_name + "_7"
        molecule_name_8 = molecule_name + "_8"
        molecule_name_9 = molecule_name + "_9"
        molecule_name_10 = molecule_name + "_10"
        molecule_name_11 = molecule_name + "_11"
        molecule_name_12 = molecule_name + "_12"
        molecule_name_13 = molecule_name + "_13"
        molecule_name_14 = molecule_name + "_14"
        molecule_name_15 = molecule_name + "_15"
        molecule_name_16 = molecule_name + "_16"
        molecule_name_17 = molecule_name + "_17"
        molecule_name_18 = molecule_name + "_18"
        molecule_name_19 = molecule_name + "_19"
        molecule_name_20 = molecule_name + "_20"
        molecule_name_21 = molecule_name + "_21"
        molecule_name_22 = molecule_name + "_22"
        molecule_name_23 = molecule_name + "_23"
        molecule_name_24 = molecule_name + "_24"
        molecule_name_25 = molecule_name + "_25"
        

        translate_line_1 = "{0.0 0.0 0.0}"
        translate_line_2 = "{0.0 0.0 " + str(z_trans) + " }"
        translate_line_3 = "{0.0 0.0 " + str(-z_trans) + " }"

        translate_line_4 = "{0.0 " + str(y_trans) + " " + str(z_trans) + "}"
        translate_line_5 = "{0.0 " + str(y_trans) + " 0.0}"
        translate_line_6 = "{0.0 " + str(y_trans) + " " + str(-z_trans) + "}"

        translate_line_7 = "{0.0 " + str(-y_trans) + " " + str(z_trans) + "}"
        translate_line_8 = "{0.0 " + str(-y_trans) + " 0.0}"
        translate_line_9 = "{0.0 " + str(-y_trans) + " " + str(-z_trans) + "}"

        translate_line_10 = "{0.0 " + str(2*y_trans) + " " + str(2*-z_trans) + "}"
        translate_line_11 = "{0.0 " + str(2*y_trans) + " " + str(-z_trans) + "}"
        translate_line_12 = "{0.0 " + str(2*y_trans) + " 0.0}"
        translate_line_13 = "{0.0 " + str(2*y_trans) + " " + str(z_trans) + "}"
        translate_line_14 = "{0.0 " + str(2*y_trans) + " " + str(2*z_trans) + "}"

        translate_line_15 = "{0.0 " + str(2*-y_trans) + " " + str(2*-z_trans) + "}"
        translate_line_16 = "{0.0 " + str(2*-y_trans) + " " + str(-z_trans) + "}"
        translate_line_17 = "{0.0 " + str(2*-y_trans) + " 0.0}"
        translate_line_18 = "{0.0 " + str(2*-y_trans) + " " + str(z_trans) + "}"
        translate_line_19 = "{0.0 " + str(2*-y_trans) + " " + str(2*z_trans) + "}"

        translate_line_20 = "{0.0 " + str(y_trans) + " " + str(2*-z_trans) + "}"
        translate_line_21 = "{0.0 " + str(y_trans) + " " + str(2*z_trans) + "}"

        translate_line_22 = "{0.0 " + str(-y_trans) + " " + str(2*-z_trans) + "}"
        translate_line_23 = "{0.0 " + str(-y_trans) + " " + str(2*z_trans) + "}"

        translate_line_24 = "{0.0 0.0 " + str(2*-z_trans) + "}"
        translate_line_25 = "{0.0 0.0 " + str(2*z_trans) + "}"

        combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + " " + molecule_name_5 + " " + molecule_name_6 + " " + molecule_name_7 + " " + molecule_name_8 + " " + molecule_name_9 + " " + molecule_name_10 + " " + molecule_name_11 + " " + molecule_name_12 + " " + molecule_name_13 + " " + molecule_name_14 + " " + molecule_name_15 + " " + molecule_name_16 + " " + molecule_name_17 + " " + molecule_name_18 + " " + molecule_name_19 + " " + molecule_name_20 + " " + molecule_name_21 + " " + molecule_name_22 + " " + molecule_name_23 + " " + molecule_name_24 + " " + molecule_name_25 + "}"

        base_mol_name = molecule_name.split("_")[0]
        intleap_path = system_name + ".intleap"

        prmtop_filepath =  os.path.join(output_dir, system_name + ".prmtop")
        rst_filepath = os.path.join(output_dir, system_name + ".rst7")
     
        three_three_array_pdb_filepath = os.path.join(output_dir, system_name + ".pdb")
 
        head_rescode, mainchain_rescode, tail_rescode = self.manager.retrieve_polymeric_rescodes(base_molecule_name)

        file_content = f"""source leaprc.gaff
        source leaprc.water.fb3
        source leaprc.protein.ff14SB

        loadamberprep {head_prepi_filepath}
        loadamberprep {mainchain_prepi_filepath}
        loadamberprep {tail_prepi_filepath}

        list
            
        {molecule_name_1} = loadpdb {pdb_file}
        {molecule_name_2} = loadpdb {pdb_file}
        {molecule_name_3} = loadpdb {pdb_file}
        {molecule_name_4} = loadpdb {pdb_file}
        {molecule_name_5} = loadpdb {pdb_file}
        {molecule_name_6} = loadpdb {pdb_file}
        {molecule_name_7} = loadpdb {pdb_file}
        {molecule_name_8} = loadpdb {pdb_file}
        {molecule_name_9} = loadpdb {pdb_file}
        {molecule_name_10} = loadpdb {pdb_file}
        {molecule_name_11} = loadpdb {pdb_file}
        {molecule_name_12} = loadpdb {pdb_file}
        {molecule_name_13} = loadpdb {pdb_file}
        {molecule_name_14} = loadpdb {pdb_file}
        {molecule_name_15} = loadpdb {pdb_file}
        {molecule_name_16} = loadpdb {pdb_file}
        {molecule_name_17} = loadpdb {pdb_file}
        {molecule_name_18} = loadpdb {pdb_file}
        {molecule_name_19} = loadpdb {pdb_file}
        {molecule_name_20} = loadpdb {pdb_file}
        {molecule_name_21} = loadpdb {pdb_file}
        {molecule_name_22} = loadpdb {pdb_file}
        {molecule_name_23} = loadpdb {pdb_file}
        {molecule_name_24} = loadpdb {pdb_file}
        {molecule_name_25} = loadpdb {pdb_file}

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
        translate {molecule_name_10} {translate_line_10}
        translate {molecule_name_11} {translate_line_11}
        translate {molecule_name_12} {translate_line_12}
        translate {molecule_name_13} {translate_line_13}
        translate {molecule_name_14} {translate_line_14}
        translate {molecule_name_15} {translate_line_15}
        translate {molecule_name_16} {translate_line_16}
        translate {molecule_name_17} {translate_line_17}
        translate {molecule_name_18} {translate_line_18}
        translate {molecule_name_19} {translate_line_19}
        translate {molecule_name_20} {translate_line_20}
        translate {molecule_name_21} {translate_line_21}
        translate {molecule_name_22} {translate_line_22}
        translate {molecule_name_23} {translate_line_23}
        translate {molecule_name_24} {translate_line_24}
        translate {molecule_name_25} {translate_line_25}
         
        system = combine {combine_line}
        setBox system vdw 0.0
        saveamberparm system {prmtop_filepath} {rst_filepath}
        savepdb system {three_three_array_pdb_filepath}
    
        quit
        """
    
        with open(intleap_path, 'w') as file:
            file.write(file_content)
            
        leap_command = "tleap -f " + intleap_path
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

        cd_command = "cd " + str(self.manager.main_dir)
        result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return(system_name)
    
    def generate_5_5_polymer_array_random(self, base_molecule_name, molecule_name):
        """
        Generates a 5x5 array of polymers using Packmol and adds 'TER' records after each polymer chain.

        This function creates a 5x5 array of polymer molecules based on a given input molecule. It uses Packmol 
        to arrange the molecules in the specified array and adds 'TER' records after each polymer chain to ensure 
        proper formatting for subsequent simulations.

        Parameters:
        dirs (object): An object that provides methods for directory and file management.
        molecule_name (str): The name of the molecule to be used for creating the polymer array. i.e. "3HB_10_polymer"
        base_molecule_name (str): The base name of the molecule used to identify the polymeric residue codes. i.e. "3HB_trimer"

        Returns:
        None

        The function performs the following steps:
        1. Defines the file subtype and constructs the PDB file path.
        2. Creates an output directory if it doesn't exist.
        3. Constructs the output PDB file name for the unsolved array.
        4. Calculates the dimensions of the box required to contain the 3x3 array of polymers.
        5. Generates the Packmol input file content with the specified parameters.
        6. Writes the Packmol input file to the output directory.
        7. Executes the Packmol command to generate the unsolved array PDB file.
        8. Adds 'TER' records to the generated PDB file using the `add_ter_to_pckml_result` function.

        Example usage:
            generate_3_3_polymer_array_pckml(dirs, 'polymer_molecule', 'base_polymer')
        """
        file_subtype = "_5_5_array_random"
        pdb_file = self.manager.load_pdb_filepath(molecule_name)
        system_name = molecule_name + file_subtype
        output_dir = os.path.join(self.manager.systems_dir, system_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        array_pdb_name = system_name + ".pdb"
        array_pdb = os.path.join(output_dir, array_pdb_name)
        x,y,z = self.get_xyz_dists(pdb_file)
        box_edge_len = int((max(x,y,z))*2)
        box_sizes = box_edge_len, box_edge_len, box_edge_len
        file_content=f"""tolerance 2.0
            output {array_pdb}
            filetype pdb
            structure {pdb_file}
            	number 25
                inside box 0. 0. 0. {box_sizes[0]}. {box_sizes[1]}. {box_sizes[2]}.
            end structure
            """
        packmol_input_file = system_name + ".pckml"
        packmol_input_filepath = os.path.join(output_dir, packmol_input_file)
        with open(packmol_input_filepath, 'w') as file:
            file.write(file_content)
        packmol_command = str(self.packmol_path) + " < " + packmol_input_filepath
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
        
        self.add_ter_to_pckml_result(array_pdb, base_molecule_name)
        return(system_name)
        
    def generate_polymer_5_5_array(self, base_molecule_name, molecule_name, method, crystal_trans=None):
        if method == "crystal":
            system_name = self.generate_5_5_polymer_array_crystal(base_molecule_name, molecule_name, crystal_trans)
        if method == "random":
            system_name = self.generate_5_5_polymer_array_random(base_molecule_name, molecule_name)
            self.gen_amber_params_4_pckml_array(system_name, base_molecule_name)
        return(system_name)
            
    def build_2_10_polymer_array(self, base_molecule_name=None, molecule_name=None, buffer=None):
         # Thsis function builds arrays of polymers using the pre generated pdb files
         if molecule_name == None or base_molecule_name == None:
            print("Please provide 2 arguments as follows: build_2_10_polymer_array(base_molecule_name, molecule_name, finite_polymer)")
            print("Base polymer name: A string of the polymer name, i.e. '3HB_trimer'")
            print("Polymer name: A string of the polymer name, i.e. '3HB_10_polymer'")
            
            return(None)        

         if buffer == None:
             buffer = "10"
         else:
             buffer = str(buffer)
        
         
         pdb_file = self.manager.load_pdb_filepath(molecule_name)
         base_pdb_file = self.manager.load_pdb_filepath(base_molecule_name)
         
         molecule_dir = os.path.join(self.manager.molecules_dir, base_molecule_name)
         cd_command = "cd " + molecule_dir
         print(cd_command)

         try:
             os.chdir(molecule_dir)
             print("Current directory:", os.getcwd())
         except Exception as e:
             print("Exception:", e)
             
         head_prepi_filepath = "head_" + base_molecule_name + ".prepi"
         mainchain_prepi_filepath = "mainchain_" + base_molecule_name + ".prepi"
         tail_prepi_filepath = "tail_" + base_molecule_name + ".prepi"

         file_subtype = "_2_10_array" 
         output_dir = os.path.join(self.manager.systems_dir, (molecule_name + file_subtype))
         if not os.path.exists(output_dir):
             os.makedirs(output_dir)

         x, y, z = self.get_xyz_dists(pdb_file)

         #translate_distance = int(max(x, y)) + 1 # Removed z as they should not overlap in this distance
         translate_distance = int((max(x, y))/2) # Removed z as they should not overlap in this distance

         molecule_name_1 = molecule_name + "_1"
         molecule_name_2 = molecule_name + "_2"
         molecule_name_3 = molecule_name + "_3"
         molecule_name_4 = molecule_name + "_4"
         molecule_name_5 = molecule_name + "_5"
         molecule_name_6 = molecule_name + "_6"
         molecule_name_7 = molecule_name + "_7"
         molecule_name_8 = molecule_name + "_8"
         molecule_name_9 = molecule_name + "_9"
         molecule_name_10 = molecule_name + "_10"
         molecule_name_11 = molecule_name + "_11"
         molecule_name_12 = molecule_name + "_12"
         molecule_name_13 = molecule_name + "_13"
         molecule_name_14 = molecule_name + "_14"
         molecule_name_15 = molecule_name + "_15"
         molecule_name_16 = molecule_name + "_16"
         molecule_name_17 = molecule_name + "_17"
         molecule_name_18 = molecule_name + "_18"
         molecule_name_19 = molecule_name + "_19"
         molecule_name_20 = molecule_name + "_20"
         
         # central 2 polymers on the top line (1,2)
         translate_line_1 = "{0.0 " + str(translate_distance/2) + " 0.0}"
         translate_line_2 = "{0.0 " + str(-translate_distance/2) + " 0.0}"

         # central 2 polymers on the bottom line (3,4)
         translate_line_3 = "{0.0 " + str(translate_distance/2) + " " + str(translate_distance) + "}"
         translate_line_4 = "{0.0 " + str(-translate_distance/2) + " " + str(translate_distance) + "}"

         # Polymers either side of central 2 polymers on the top line (5,6)
         translate_line_5 = "{0.0 " + str(translate_distance*1.5) + " 0.0}"
         translate_line_6 = "{0.0 " + str(-translate_distance*1.5) + " 0.0}"

         # Polymers either side of central 2 polymers on the bottom line (7,8)
         translate_line_7 = "{0.0 " + str(translate_distance*1.5) + " " + str(translate_distance) + "}"
         translate_line_8 = "{0.0 " + str(-translate_distance*1.5) + " " + str(translate_distance) + "}"

         # Polymers (9,10)
         translate_line_9 = "{0.0 " + str(translate_distance*2.5) + " 0.0}"
         translate_line_10 = "{0.0 " + str(-translate_distance*2.5) + " 0.0}"

         # Polymers (11,12)
         translate_line_11 = "{0.0 " + str(translate_distance*2.5) + " " + str(translate_distance) + "}"
         translate_line_12 = "{0.0 " + str(-translate_distance*2.5) + " " + str(translate_distance) + "}"

         # Polymers (13,14)
         translate_line_13 = "{0.0 " + str(translate_distance*3.5) + " 0.0}"
         translate_line_14 = "{0.0 " + str(-translate_distance*3.5) + " 0.0}"

         # Polymers (15,16)
         translate_line_15 = "{0.0 " + str(translate_distance*3.5) + " " + str(translate_distance) + "}"
         translate_line_16 = "{0.0 " + str(-translate_distance*3.5) + " " + str(translate_distance) + "}"

         # Polymers (17,18)
         translate_line_17 = "{0.0 " + str(translate_distance*4.5) + " 0.0}"
         translate_line_18 = "{0.0 " + str(-translate_distance*4.5) + " 0.0}"

         # Polymers (19,20)
         translate_line_19 = "{0.0 " + str(translate_distance*4.5) + " " + str(translate_distance) + "}"
         translate_line_20 = "{0.0 " + str(-translate_distance*4.5) + " " + str(translate_distance) + "}"

         combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + " " + molecule_name_5 + " " + molecule_name_6 + " " + molecule_name_7 + " " + molecule_name_8 + " " + molecule_name_9 + " " + molecule_name_10 + " " + molecule_name_11 + " " + molecule_name_12 + " " + molecule_name_13 + " " + molecule_name_14 + " " + molecule_name_15 + " " + molecule_name_16 + " " + molecule_name_17 + " " + molecule_name_18 + " " + molecule_name_19 + " " + molecule_name_20 + "}"

         base_mol_name = molecule_name.split("_")[0]
         intleap_path = base_mol_name + file_subtype + ".intleap"

         system_name = molecule_name + file_subtype
         unsolved_system_name = "unsolved_" + molecule_name + file_subtype
        
         prmtop_filepath =  os.path.join(output_dir, system_name + ".prmtop")
         rst_filepath = os.path.join(output_dir, system_name + ".rst7")

         unsolved_prmtop_filepath =  os.path.join(output_dir, "unsolved_" + molecule_name + file_subtype + ".prmtop")
         unsolved_rst_filepath = os.path.join(output_dir, "unsolved_" + molecule_name + file_subtype + ".rst7")
        
         two_ten_array_pdb_filepath = os.path.join(output_dir, system_name + ".pdb")
         unsolved_two_ten_array_pdb_filepath = os.path.join(output_dir, "unsolved_" + molecule_name + file_subtype + ".pdb")
    
         head_rescode, mainchain_rescode, tail_rescode = self.manager.retrieve_polymeric_rescodes(base_molecule_name)

         file_content = f"""source leaprc.gaff
         source leaprc.water.fb3
         source leaprc.protein.ff14SB

         loadamberprep {head_prepi_filepath}
         loadamberprep {mainchain_prepi_filepath}
         loadamberprep {tail_prepi_filepath}

         list
            
         {molecule_name_1} = loadpdb {pdb_file}
         {molecule_name_2} = loadpdb {pdb_file}
         {molecule_name_3} = loadpdb {pdb_file}
         {molecule_name_4} = loadpdb {pdb_file}
         {molecule_name_5} = loadpdb {pdb_file}
         {molecule_name_6} = loadpdb {pdb_file}
         {molecule_name_7} = loadpdb {pdb_file}
         {molecule_name_8} = loadpdb {pdb_file}
         {molecule_name_9} = loadpdb {pdb_file}
         {molecule_name_10} = loadpdb {pdb_file}
         {molecule_name_11} = loadpdb {pdb_file}
         {molecule_name_12} = loadpdb {pdb_file}
         {molecule_name_13} = loadpdb {pdb_file}
         {molecule_name_14} = loadpdb {pdb_file}
         {molecule_name_15} = loadpdb {pdb_file}
         {molecule_name_16} = loadpdb {pdb_file}
         {molecule_name_17} = loadpdb {pdb_file}
         {molecule_name_18} = loadpdb {pdb_file}
         {molecule_name_19} = loadpdb {pdb_file}
         {molecule_name_20} = loadpdb {pdb_file}
         
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
         translate {molecule_name_10} {translate_line_10}
         translate {molecule_name_11} {translate_line_11}
         translate {molecule_name_12} {translate_line_12}
         translate {molecule_name_13} {translate_line_13}
         translate {molecule_name_14} {translate_line_14}
         translate {molecule_name_15} {translate_line_15}
         translate {molecule_name_16} {translate_line_16}
         translate {molecule_name_17} {translate_line_17}
         translate {molecule_name_18} {translate_line_18}
         translate {molecule_name_19} {translate_line_19}
         translate {molecule_name_20} {translate_line_20}

         system = combine {combine_line}
         unsolved_system = system
         setBox unsolved_system vdw 0.0
         saveamberparm unsolved_system {unsolved_prmtop_filepath} {unsolved_rst_filepath}
         savepdb unsolved_system {unsolved_two_ten_array_pdb_filepath}
    
         solvatebox system TIP3PBOX {buffer}

         saveamberparm system {prmtop_filepath} {rst_filepath}
         savepdb system {two_ten_array_pdb_filepath}
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

         cd_command = "cd " + str(self.manager.main_dir)
         result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
         return(system_name, unsolved_system_name)

    def build_2_10_polymer_array_infinite(self, base_molecule_name=None, molecule_name=None, number_of_units=None, buffer=None):
        # Thsis function builds arrays of polymers using the pre generated pdb files
        if number_of_units == None or base_molecule_name == None or molecule_name == None:
            print("Please provide 3 arguments as follows: build_2_10_polymer_array(base_molecule_name, polymer_name, number_of_units)")
            print("Base polymer name: A string of the polymer name, i.e. '3HB_trimer'")
            print("Polymer name: A string of the polymer name, i.e. '3HB_10_polymer'")
            
            return(None)        

        if buffer == None:
            buffer = None
        else:
            buffer = str(buffer)
        
         
        pdb_file = self.manager.load_pdb_filepath(molecule_name)
         
        molecule_dir = os.path.join(self.manager.molecules_dir, base_molecule_name)
        cd_command = "cd " + molecule_dir
        print(cd_command)

        try:
            os.chdir(molecule_dir)
            print("Current directory:", os.getcwd())
        except Exception as e:
            print("Exception:", e)
             
        mainchain_prepi_filepath = "mainchain_" + base_molecule_name + ".prepi"

        file_subtype = "_2_10_array_infinite" 
        output_dir = os.path.join(self.manager.systems_dir, (molecule_name + file_subtype))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        x, y, z = self.get_xyz_dists(pdb_file)

        #translate_distance = int(max(x, y)) + 1 # Removed z as they should not overlap in this distance
        translate_distance = int((max(x, y))/2) # Removed z as they should not overlap in this distance

        molecule_name_1 = molecule_name + "_1"
        molecule_name_2 = molecule_name + "_2"
        molecule_name_3 = molecule_name + "_3"
        molecule_name_4 = molecule_name + "_4"
        molecule_name_5 = molecule_name + "_5"
        molecule_name_6 = molecule_name + "_6"
        molecule_name_7 = molecule_name + "_7"
        molecule_name_8 = molecule_name + "_8"
        molecule_name_9 = molecule_name + "_9"
        molecule_name_10 = molecule_name + "_10"
        molecule_name_11 = molecule_name + "_11"
        molecule_name_12 = molecule_name + "_12"
        molecule_name_13 = molecule_name + "_13"
        molecule_name_14 = molecule_name + "_14"
        molecule_name_15 = molecule_name + "_15"
        molecule_name_16 = molecule_name + "_16"
        molecule_name_17 = molecule_name + "_17"
        molecule_name_18 = molecule_name + "_18"
        molecule_name_19 = molecule_name + "_19"
        molecule_name_20 = molecule_name + "_20"
         
        # central 2 polymers on the top line (1,2)
        translate_line_1 = "{0.0 " + str(translate_distance/2) + " 0.0}"
        translate_line_2 = "{0.0 " + str(-translate_distance/2) + " 0.0}"

        # central 2 polymers on the bottom line (3,4)
        translate_line_3 = "{0.0 " + str(translate_distance/2) + " " + str(translate_distance) + "}"
        translate_line_4 = "{0.0 " + str(-translate_distance/2) + " " + str(translate_distance) + "}"

        # Polymers either side of central 2 polymers on the top line (5,6)
        translate_line_5 = "{0.0 " + str(translate_distance*1.5) + " 0.0}"
        translate_line_6 = "{0.0 " + str(-translate_distance*1.5) + " 0.0}"

        # Polymers either side of central 2 polymers on the bottom line (7,8)
        translate_line_7 = "{0.0 " + str(translate_distance*1.5) + " " + str(translate_distance) + "}"
        translate_line_8 = "{0.0 " + str(-translate_distance*1.5) + " " + str(translate_distance) + "}"

        # Polymers (9,10)
        translate_line_9 = "{0.0 " + str(translate_distance*2.5) + " 0.0}"
        translate_line_10 = "{0.0 " + str(-translate_distance*2.5) + " 0.0}"

        # Polymers (11,12)
        translate_line_11 = "{0.0 " + str(translate_distance*2.5) + " " + str(translate_distance) + "}"
        translate_line_12 = "{0.0 " + str(-translate_distance*2.5) + " " + str(translate_distance) + "}"

        # Polymers (13,14)
        translate_line_13 = "{0.0 " + str(translate_distance*3.5) + " 0.0}"
        translate_line_14 = "{0.0 " + str(-translate_distance*3.5) + " 0.0}"

        # Polymers (15,16)
        translate_line_15 = "{0.0 " + str(translate_distance*3.5) + " " + str(translate_distance) + "}"
        translate_line_16 = "{0.0 " + str(-translate_distance*3.5) + " " + str(translate_distance) + "}"

        # Polymers (17,18)
        translate_line_17 = "{0.0 " + str(translate_distance*4.5) + " 0.0}"
        translate_line_18 = "{0.0 " + str(-translate_distance*4.5) + " 0.0}"

        # Polymers (19,20)
        translate_line_19 = "{0.0 " + str(translate_distance*4.5) + " " + str(translate_distance) + "}"
        translate_line_20 = "{0.0 " + str(-translate_distance*4.5) + " " + str(translate_distance) + "}"

        combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + " " + molecule_name_5 + " " + molecule_name_6 + " " + molecule_name_7 + " " + molecule_name_8 + " " + molecule_name_9 + " " + molecule_name_10 + " " + molecule_name_11 + " " + molecule_name_12 + " " + molecule_name_13 + " " + molecule_name_14 + " " + molecule_name_15 + " " + molecule_name_16 + " " + molecule_name_17 + " " + molecule_name_18 + " " + molecule_name_19 + " " + molecule_name_20 + "}"

        base_mol_name = molecule_name.split("_")[0]
        intleap_path = base_mol_name + file_subtype + ".intleap"

        system_name = molecule_name + file_subtype
        unsolved_system_name = "unsolved_" + molecule_name + file_subtype
        
        prmtop_filepath =  os.path.join(output_dir, system_name + ".prmtop")
        rst_filepath = os.path.join(output_dir, system_name + ".rst7")

        unsolved_prmtop_filepath =  os.path.join(output_dir, "unsolved_" + molecule_name + file_subtype + ".prmtop")
        unsolved_rst_filepath = os.path.join(output_dir, "unsolved_" + molecule_name + file_subtype + ".rst7")
        
        two_ten_array_pdb_filepath = os.path.join(output_dir, system_name + ".pdb")
        unsolved_two_ten_array_pdb_filepath = os.path.join(output_dir, "unsolved_" + molecule_name + file_subtype + ".pdb")
    
        mainchain_rescode = self.manager.retrieve_polymeric_rescodes(base_molecule_name)[1]

        poly_string = " ".join([mainchain_rescode] * number_of_units)
        poly_string = f"{{{poly_string}}}"

        file_content = f"""source leaprc.gaff
        source leaprc.water.fb3
        source leaprc.protein.ff14SB

        loadamberprep {mainchain_prepi_filepath}

        list
            
        {molecule_name_1} = sequence {poly_string}
        {molecule_name_2} = sequence {poly_string}
        {molecule_name_3} = sequence {poly_string}
        {molecule_name_4} = sequence {poly_string}
        {molecule_name_5} = sequence {poly_string}
        {molecule_name_6} = sequence {poly_string}
        {molecule_name_7} = sequence {poly_string}
        {molecule_name_8} = sequence {poly_string}
        {molecule_name_9} = sequence {poly_string}
        {molecule_name_10} = sequence {poly_string}
        {molecule_name_11} = sequence {poly_string}
        {molecule_name_12} = sequence {poly_string}
        {molecule_name_13} = sequence {poly_string}
        {molecule_name_14} = sequence {poly_string}
        {molecule_name_15} = sequence {poly_string}
        {molecule_name_16} = sequence {poly_string}
        {molecule_name_17} = sequence {poly_string}
        {molecule_name_18} = sequence {poly_string}
        {molecule_name_19} = sequence {poly_string}
        {molecule_name_20} = sequence {poly_string}
         
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
        translate {molecule_name_10} {translate_line_10}
        translate {molecule_name_11} {translate_line_11}
        translate {molecule_name_12} {translate_line_12}
        translate {molecule_name_13} {translate_line_13}
        translate {molecule_name_14} {translate_line_14}
        translate {molecule_name_15} {translate_line_15}
        translate {molecule_name_16} {translate_line_16}
        translate {molecule_name_17} {translate_line_17}
        translate {molecule_name_18} {translate_line_18}
        translate {molecule_name_19} {translate_line_19}
        translate {molecule_name_20} {translate_line_20}

        system = combine {combine_line}
        unsolved_system = system
        setBox unsolved_system vdw 0.0
        saveamberparm unsolved_system {unsolved_prmtop_filepath} {unsolved_rst_filepath}
        savepdb unsolved_system {unsolved_two_ten_array_pdb_filepath}

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

        cd_command = "cd " + str(self.manager.main_dir)
        result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return(system_name, unsolved_system_name)

    def solvate_frame_from_sim(self, base_molecule_name, system_name, pdb_file, box):
    
        molecule_dir = os.path.join(self.manager.molecules_dir, base_molecule_name)
        cd_command = "cd " + molecule_dir
        try:
            os.chdir(molecule_dir)
        except Exception as e:
            print("Exception:", e)
    
        mainchain_prepi_filepath = "mainchain_" + base_molecule_name + ".prepi"
        head_prepi_filepath = "head_" + base_molecule_name + ".prepi"
        tail_prepi_filepath = "tail_" + base_molecule_name + ".prepi"
        
        output_dir = os.path.join(self.manager.systems_dir, system_name)
        print(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        x_vec = box[0]
        y_vec = box[1]
        z_vec = box[2]

        box_set_line = "set system box {" + str(int(x_vec)+1) + " " + str(int(y_vec)+1) + " " + str(int(z_vec)+1) + " 90.0 90.0 90.0}"

        intleap_path = system_name + ".intleap"
        prmtop_path = os.path.join(output_dir, system_name + ".prmtop")
        rst_path = os.path.join(output_dir, system_name + "rst7")
        solvated_pdb = os.path.join(output_dir, system_name + ".pdb")
    
    
        file_content = f"""source leaprc.gaff
        source leaprc.water.fb3
    
        loadamberprep {head_prepi_filepath}
        loadamberprep {mainchain_prepi_filepath}
        loadamberprep {tail_prepi_filepath}
        list
        system = loadpdb {pdb_file}
        solvatebox system TIP3PBOX 0.0

        saveamberparm system {prmtop_path} {rst_path}
        savepdb system {solvated_pdb}
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

        cd_command = "cd " + str(self.manager.main_dir)
        result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return(system_name)
    
    def gen_amber_params_sing_mol(self, molecule_name, box_radius=None):
        """
        Generates a polymer PDB file using `tleap` based on the specified molecule and number of units.

        Args:
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
        if box_radius==None:
            box_radius = 10.0
        else:
            box_radius = box_radius
        
        molecule_prepi_file = os.path.join(self.manager.molecules_dir, molecule_name, (molecule_name + ".prepi"))
        molecule_frcmod_file = os.path.join(self.manager.molecules_dir, molecule_name, (molecule_name + ".frcmod"))

        file_subtype = "_sing_mol"
        output_name = molecule_name + file_subtype

        output_dir = os.path.join(self.manager.systems_dir, output_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)


        intleap_path = os.path.join(output_dir, output_name + ".intleap")
        prmtop_filepath = os.path.join(output_dir, output_name + ".prmtop")
        rst_filepath = os.path.join(output_dir, output_name + ".rst7")
        pdb_filepath = os.path.join(output_dir, output_name + ".pdb")

        rescode = self.manager.retrieve_rescode(molecule_name)

        file_content = f"""source leaprc.gaff
             source leaprc.water.fb3
             source leaprc.protein.ff14SB

             loadamberprep {molecule_prepi_file}
             loadamberparams {molecule_frcmod_file}

             list

             system = {rescode}
             setBox system centers {box_radius}
             saveamberparm system {prmtop_filepath} {rst_filepath}
             savepdb system {pdb_filepath}
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
        return(output_name)  

        
#class BuildBioOilSystems(BuildSystems):
#    def __init__:
 #       pass

    