# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 11:46:42 2024

@author: danie
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolFromPDBFile
import subprocess
import numpy as np
import sys
import os
'''
This script uses the files generated with "python_parameterizer.py" to construct a tleap file for 
    a system where an array of 4 molecules are solvated.  

The size of these systems is the maximum distance between any 2 atoms in the molecule being
    solvated multiplied by 4 (as we have more molecules in the system):
        
        box_dimensions = (x, y, z) ; where x, y and z == (max_pair_distance(mol))*2

1. Get the current working directory as a variable called "current_working_directory"
    Note: this uses bash commands that are executed with python - 
        run the command in the unbuntu or SCW shell to see what it does.
'''
# Bash command to print the current working directory (pwd)
pwd_command = "pwd"
# Run the command and capture the output
result = subprocess.run(pwd_command, shell=True, stdout=subprocess.PIPE, text=True)
# Extract the output into a variable
current_working_directory = result.stdout.strip()
'''
2. Construct a filepath to the molecule directory:
   i.e. ~/polymer_md_1/pdb_files/molecules/ethane
   The molecule name is provided as a command line argument and the path is constructed to it.
   If the directroy doesnt exist the script will exit and ask you to paramerterize your molecule.
'''
molecule_name = sys.argv[1]
molecule_path = current_working_directory + "/pdb_files/molecules"
pdb_directory_path = os.path.join(molecule_path, molecule_name)

# Check if the directory exists
if not os.path.exists(pdb_directory_path):
    print(f"Error: The directory '{pdb_directory_path}' does not exist.")
    print(f"Please parameterize '{molecule_name}' from its pdb file using 'python_paramertizer.py <path/to/{molecule_name}.pdb>'")
    sys.exit(1)  # Exit with a non-zero code to indicate an error
'''
3. Construct filepaths to the files we need to create the solvated molecule
   Construct filepaths we will be making (i.e. .rst7 and .prmtop)
'''
pdb_filepath = os.path.join(pdb_directory_path, (molecule_name + ".pdb"))
mol2_filepath = os.path.join(pdb_directory_path, (molecule_name + ".mol2"))
prmtop_filepath = os.path.join(pdb_directory_path, (molecule_name + "_2_2_array.prmtop"))
rst_filepath = os.path.join(pdb_directory_path, (molecule_name + "_2_2_array.rst7"))
two_two_array_pdb_filepath = os.path.join(pdb_directory_path, (molecule_name + "_2_2_array.pdb"))
'''
4. Load the molecule from the pdb file using rdkit
'''
mol = MolFromPDBFile(pdb_filepath)
'''
5. Define function that finds the maximum distance between pairs of atoms in the molecule
'''
def max_pairwise_distance(mol):
    """
    Calculate the maximum pairwise distance between atoms in a molecule.

    Parameters:
    - mol (Chem.Mol): RDKit molecule object.

    Returns:
    - float: Maximum pairwise distance between atoms.

    This function extracts the atomic coordinates from the provided RDKit molecule object,
    calculates the pairwise distances between atoms, and returns the maximum distance.
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
'''
6. Calculate max distance between 2 parts of a molecule and specify the box distance as twice this
'''
max_dist = max_pairwise_distance(mol)

# We will translate the molecules from the centre by their max pariwise distance.
# This will give an overall distance of (2*max_pairwise_distance) between each molecule
translate_dist = float((int(max_dist)+1))

# The box_distance is (4*max_pairwise_distance) as we have 4 molecules
box_dist = float((int(max_dist)+1)*4) # the +1 is here as the "int" function converts a decimal to an integer but always rounds down
'''
7. Generate an intleap file as seen in the jupyter notebook "python_parameterizer.py"

    This section uses filepaths defined in step 3.

Note: You can edit the contents of the file, but I recommend testing the commands
    you wish to execute by loading tleap in the shell and entering each line manually.
'''
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

# Filepath for the intleap file
intleap_path = os.path.join(pdb_directory_path, (molecule_name + ".2_2_array_intleap"))

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
'''
8. Execute the tleap script with python subprocess module
''' 
leap_command = "tleap -f " + intleap_path
subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


