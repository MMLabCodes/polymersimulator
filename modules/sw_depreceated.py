# This file contains functions that are no longer used but can still be called as legacy functions.
# Support will not be continued for these functions.

def build_3_3_polymer_array_old(self, directories=None, molecule_name=None, number_of_units=None):
    # This is an old function that builds 3_3_arrays of polymers using tleap
    # Each polymer is generated in the tleap script
    # Issues with the function:
    #    - Final box too big
    #    - Some polymers not are not in the correct conformation and too large forces mean simulations do not run
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

    file_subtype = "_3_3_array_" + str(number_of_units) + "_polymer"
    
    head_prepi_filepath = "head_" + molecule_name + ".prepi"
    mainchain_prepi_filepath = "mainchain_" + molecule_name + ".prepi" 
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
    translate_line_3 = "{" + str(-translate_distance) + " 0.0 0.0}"
        

    translate_line_4 = "{" + str(translate_distance) + " " + str(translate_distance) + " 0.0}"       
    translate_line_5 = "{" + str(-translate_distance) + " " + str(translate_distance) + " 0.0}" 
    translate_line_6 = "{0.0 " + str(translate_distance) + " 0.0}"

    translate_line_7 = "{" + str(translate_distance) + " " + str(-translate_distance) + " 0.0}"      
    translate_line_8 = "{" + str(-translate_distance) + " " + str(-translate_distance) + " 0.0}"      
    translate_line_9 = "{0.0 " + str(-translate_distance) + " 0.0}"

    combine_line = "{" + molecule_name_1 + " " + molecule_name_2 + " " + molecule_name_3 + " " + molecule_name_4 + " " + molecule_name_5 + " " + molecule_name_6 + " " + molecule_name_7 + " " + molecule_name_8 + " " + molecule_name_9 + "}"

    base_mol_name = molecule_name.split("_")[0]     
    intleap_path = base_mol_name + file_subtype + ".intleap"
    
    prmtop_filepath =  os.path.join(output_dir, base_mol_name + file_subtype + "_" + str(box_dist) + ".prmtop")
    rst_filepath = os.path.join(output_dir, base_mol_name + file_subtype + "_" + str(box_dist) + ".rst7")

    unsolved_prmtop_filepath =  os.path.join(output_dir, "unsolved_" + base_mol_name + file_subtype + ".prmtop")
    unsolved_rst_filepath = os.path.join(output_dir, "unsolved_" + base_mol_name + file_subtype + ".rst7")
        
    three_three_array_pdb_filepath = os.path.join(output_dir, base_mol_name + file_subtype + "_" + str(box_dist) + ".pdb")
    unsolved_three_three_array_pdb_filepath = os.path.join(output_dir, "unsolved_" + base_mol_name + file_subtype + ".pdb")
    
    head_rescode, mainchain_rescode, tail_rescode = directories.retrieve_polymeric_rescodes("3HB_trimer")
    
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
    print(intleap_path)
    try:
        result = subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:            
            print("Output:", result.stdout)
        else:          
            print("Error:", result.stderr)
    except Exception as e:
        # Exception occurred during subprocess execution
        print("Exception:", e)

    cd_command = "cd " + str(directories.main_dir)
    result = subprocess.run(cd_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return()

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
    subprocess.run(leap_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

def parametrized_mols_avail(self):
    '''
    Old function to retrieve molecules with mol2 files.
    Replace by the method "SnippetSimManage.mol2_avail()"
    '''
    a = False
    for root, dirs, files in os.walk(self.molecules_dir):
        dirs[:] = [d for d in dirs if d != 'depreceated']
        # Check each file in the current directory
        for file in files:
            # Check if the file has a .pdb extension
            if file.endswith(".mol2"):
                a = True
                # Construct the full path to the .pdb file
                pdb_file_path = os.path.join(root, file)
                # Extract molecule name
                pdb_file = pdb_file_path.split("/")[-1]
                print(pdb_file)
    if a == False:
        print("No parametrized molecules.")

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

def solvate_polymer_pdb(self, molecule_name, polymer_name, buffer=None):
    # This function will solvate a single polymer
    # Define directory where prepi files are found
    molecule_dir = os.path.join(self.manager.molecules_dir, molecule_name)

    if buffer == None:
         buffer = "10"
    else:
         buffer = str(buffer)

     # Change to prepi file directory
    cd_command = "cd " + molecule_dir
    try:
        os.chdir(molecule_dir)
        print("Current directory:", os.getcwd())
    except Exception as e:
           print("Exception:", e)

    # Define prepi filepaths
    head_prepi_filepath = "head_" + molecule_name + ".prepi"
    mainchain_prepi_filepath = "mainchain_" + molecule_name + ".prepi"
    tail_prepi_filepath = "tail_" + molecule_name + ".prepi"

    file_subtype = "_wat_solv"
    # Define output directory - based on polymer name
    output_dir = os.path.join(self.manager.systems_dir, polymer_name + file_subtype)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Redundant lol
    mol_name_1 = polymer_name

    # Get box size by taking x coord (its the length of the polymer) and adding a buffer
    pdb_filepath = os.path.join(self.manager.systems_dir, polymer_name, polymer_name + ".pdb")

    # Define required filepaths   
    intleap_path = polymer_name + file_subtype + ".intleap"
    filename = polymer_name + file_subtype + f"_{buffer}"
    prmtop_filepath =  os.path.join(output_dir, filename + ".prmtop")
    rst_filepath = os.path.join(output_dir, filename + ".rst7")
    pdb_outpath = os.path.join(output_dir, filename + ".pdb") 
    single_chain_pdb_filepath = os.path.join(output_dir, filename + ".pdb")


    # Define file content for leap file
    file_content = f"""source leaprc.gaff
            source leaprc.water.fb3
            source leaprc.protein.ff14SB

            loadamberprep {head_prepi_filepath}
            loadamberprep {mainchain_prepi_filepath}
            loadamberprep {tail_prepi_filepath}

            {mol_name_1} = loadpdb {pdb_filepath}

            check {mol_name_1}

            solvatebox {mol_name_1} TIP3PBOX {buffer}

            saveamberparm {mol_name_1} {prmtop_filepath} {rst_filepath}
            savepdb {mol_name_1} {pdb_outpath}
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
    return(filename)
    


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


    def add_ter_to_pckml_result(self, pdb_file, base_molecule_name):
        """
        Adds 'TER' records to a PDB file generated by Packmol after each polymer chain.

        This function is intended to be used after Packmol has built systems of polymers with different
        head, main, and tail units. Packmol does not add 'TER' records after each residue, making it 
        impossible to construct a parameter file. This function remedies that by adding 'TER' after each 
        polymer chain in the Packmol file.

        Parameters:
        dirs (object): An object that provides a method `retrieve_polymeric_rescodes` to retrieve 
                   residue codes for head, mainchain, and tail units.
        pdb_file (str): The path to the PDB file generated by Packmol.
        base_molecule_name (str): The base name of the molecule used to identify the polymeric residue codes.

        Returns:
        None

        The function performs the following steps:
        1. Retrieves the residue codes for the head, mainchain, and tail units.
        2. Reads the contents of the specified PDB file.
        3. Iterates through each line in the file to identify points where 'TER' records should be added:
            - Adds a 'TER' record when transitioning from a tail residue to a head residue.
            - Adds a 'TER' record before the 'END' record in the PDB file.
        4. Writes the modified contents back to the PDB file with the appropriate 'TER' records included.

        Example usage:
            add_ter_to_pckml_result(dirs, 'polymer_system.pdb', '3HB_trimer')

        Note:
            A situation where this function is used by itself is not present and it is used as part of the 'generate_3_3_polymer_array_pckml'
        """
        # Retrieve residue codes
        head_rescode, mainchain_rescode, tail_rescode = self.manager.retrieve_polymeric_rescodes(base_molecule_name)

        # Open file
        with open(pdb_file) as file:
            lines = file.readlines()

        # Initiate an empty list foe new lines and set the residue code to None before itereating
        modified_lines = []
        previous_residue_code = None

        # Iterate over each line in the file
        for line in lines:
            # Split the line into columns based on spaces
            columns = line.split()

            if len(columns) > 3:
                # Extract the current residue code (4th column in this case)
                current_residue_code = columns[3]

                # Check if the previous residue code was 'tAD' and the current is 'hAD'
                if previous_residue_code == tail_rescode and current_residue_code == head_rescode:
                    # Append the "TER" line to the modified lines list
                    modified_lines.append("TER\n")

                # Update the previous residue code
                previous_residue_code = current_residue_code

            if "END" in columns:
                modified_lines.append("TER\n")
    
            # Append the current line to the modified lines list
            modified_lines.append(line)

        with open(pdb_file, 'w') as file:
            file.writelines(modified_lines)
        return(None)

def gen_amber_params_sing_poly(self, base_molecule_name, polymer_name, box_radius=None):

    head_prepi_file = os.path.join(self.manager.molecules_dir, base_molecule_name, ("head_" + base_molecule_name + ".prepi"))
    mainchain_prepi_file = os.path.join(self.manager.molecules_dir, base_molecule_name, ("mainchain_" + base_molecule_name + ".prepi"))
    tail_prepi_file = os.path.join(self.manager.molecules_dir, base_molecule_name, ("tail_" + base_molecule_name + ".prepi"))
    base_molecule_frcmod_file = os.path.join(self.manager.molecules_dir, base_molecule_name, (base_molecule_name + ".frcmod"))
    polymer_pdb_file = os.path.join(self.manager.systems_dir, polymer_name, (polymer_name + ".pdb"))

    file_subtype = "_sing_poly"
    output_name = polymer_name + file_subtype

    output_dir = os.path.join(self.manager.systems_dir, output_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    intleap_path = os.path.join(output_dir, output_name + ".intleap")
    prmtop_filepath = os.path.join(output_dir, output_name + ".prmtop")
    rst_filepath = os.path.join(output_dir, output_name + ".rst7")
    pdb_filepath = os.path.join(output_dir, output_name + ".pdb")

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
        
    file_content = f"""source leaprc.gaff
             source leaprc.water.fb3
             source leaprc.protein.ff14SB

             loadamberprep {head_prepi_file}
             loadamberprep {mainchain_prepi_file}
             loadamberprep {tail_prepi_file}
             loadamberparams {base_molecule_frcmod_file}

             list

             system = loadpdb {polymer_pdb_file}
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
























