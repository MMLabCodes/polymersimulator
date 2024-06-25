# This file contains functions that are no longer used but can still be called as legacy functions.
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



































