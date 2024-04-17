# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 09:50:22 2024

@author: danie
"""

## IMPORTANT - this is now integrated into "sw_build_systems"

import subprocess
import os
import sys
import shutil
'''
This script uses a pdb file and generates a series of parameter files for the molecule in
    the pdb file for use in with scripts that build systems from molecule parameters.

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
# Initial copying and set up of directory
'''
2. Construct the path the pdb_file directory   
'''
pdb_file_destination = current_working_directory + "/pdb_files/molecules"
'''
3. Define funtions to remove duplicate bond info

Note: In **parametrization.ipynb** the need for this function is shown. Double bonds are interpreted as errors
    and instead the pdb file must be editted to only contain single bonds. This is not a problem though as the 
    positions are the same, and AMBER will still assign the correct atom types to double bonded groups (i.e. carbonyls)
    based on neighbouring atoms and bonds that do exist with that atom.
'''
def pdb_parameter_modified(pdb_file): # pdb_file is the entire path to the pdb_file - see jupyter notebook
    """
    Modify a PDB file by removing duplicate bond information.

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
'''
4. Run this function on the pdb file that parameters will be generated for
'''
input_pdb = sys.argv[1]
pdb_parameter_modified(input_pdb)
'''
5. Copy the pdb_file into a new directory - we make a new directory to store all related files
       to this molecule - otherwise things will get messy quick when we generate parameters for many molecules.
   Then we remove the original pdb_file as it now exists in its own directory.
'''
pdb_name = os.path.basename(input_pdb).split(".")[0]   
molecule_dir_path = os.path.join(pdb_file_destination, pdb_name)
os.makedirs(molecule_dir_path, exist_ok=True)

## Code from here - input in sw_build_systems.py

shutil.copy2(input_pdb, molecule_dir_path)
os.remove(input_pdb)
'''
5. Now we construct some filepaths for pdb and mol2 files
      And run the antechamber command using the subporcess python module to make this mol2 file.
'''
pdb_path = os.path.join(molecule_dir_path, (pdb_name + ".pdb"))
mol2_path = os.path.join(molecule_dir_path, (pdb_name + ".mol2"))
antechamber_command = "antechamber -i " + pdb_path + " -fi pdb -o " + mol2_path + " -fo mol2 -c bcc -s 2"
subprocess.run(antechamber_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
'''
6. Now we construct a filepath for a frcmod file
       And run the parmchk2 command using the previoulsy constructed mol2 file and the new frcmod filepath
'''
frcmod_path = os.path.join(molecule_dir_path, (pdb_name + ".frcmod"))
parmchk_command = "parmchk2 -i " + mol2_path + " -f mol2 -o " + frcmod_path
subprocess.run(parmchk_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
'''
7. Now we construct a filepath for the lib file and intleap file
   And write the inleap file with these constructed variables implemented
'''
lib_path = os.path.join(molecule_dir_path, (pdb_name + ".lib"))
intleap_path = os.path.join(molecule_dir_path, (pdb_name + ".intleap"))
file_content = f"""source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.fb3
{pdb_name} = loadmol2 {mol2_path}
saveoff {pdb_name} {lib_path}
quit
"""
with open(intleap_path, 'w') as file:
    file.write(file_content)
leap_command = "tleap -f " + intleap_path
'''
8. Run the intleap file with the python subprocess module
'''
subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
