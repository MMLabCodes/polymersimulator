from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import math

def vol_from_smiles(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol) 
    return volume

def vol_from_mol(mol):
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol)   
    return volume

def estimated_volume(pdb_file_path):
    # Van der Waals radii (in angstroms) for common elements
    # https://periodictable.com/Properties/A/VanDerWaalsRadius.v.html
    vdw_radii = {
        'H': 1.20,  # Hydrogen
        'C': 1.70,  # Carbon
        'N': 1.55,  # Nitrogen
        'O': 1.52,  # Oxygen
        'P': 1.80,  # Phosphorus
        'S': 1.80,  # Sulfur
        'F': 1.47,  # Fluorine
        'Cl': 1.75, # Chlorine
        'Br': 1.85, # Bromine
        'I': 1.98   # Iodine
    }

    # Load the PDB file
    mol = Chem.MolFromPDBFile(pdb_file_path, removeHs=False)
    if not mol:
        print("Failed to load molecule from PDB.")
        return None

    # Initialize total volume
    total_volume = 0.0

    # Iterate over atoms to sum up their volumes
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        
        # Check if the atom's VDW radius is in our dictionary
        if symbol in vdw_radii:
            radius = vdw_radii[symbol]
            # Calculate volume of the atom and add to total volume
            atom_volume = (4/3) * math.pi * (radius ** 3)
            total_volume += atom_volume
        else:
            print(f"Warning: Van der Waals radius not found for atom type '{symbol}'.")

    return total_volume/2

# Function to calculate volume from a PDB file
def vol_from_pdb(pdb_file):
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)  # Load the molecule from PDB file
    if mol is None:
        raise ValueError(f"Could not load molecule from PDB file: {pdb_file}")
    
    mol = Chem.AddHs(mol)  # Ensure hydrogens are added for volume calculation
    AllChem.EmbedMolecule(mol)  # Generate 3D coordinates for the molecule
    volume = AllChem.ComputeMolVolume(mol)  # Calculate the molecular volume
    return volume

def has_heteroatoms(mol):
    heteroatom_smarts = ["[#7]", "[#8]", "[#16]"]
    heteroatom_names = ["nitrogen", "oxygen", "sulfur"]
    heteroatoms_in_mol = []
    for i in range(len(heteroatom_smarts)):
                   pattern_smarts = Chem.MolFromSmarts(heteroatom_smarts[i])
                   if mol.HasSubstructMatch(pattern_smarts) == True:
                       heteroatoms_in_mol.append(heteroatom_names[i])
    if len(heteroatoms_in_mol) == 0:
        heteroatoms_in_mol.append("No heteroatoms")
    return(heteroatoms_in_mol)

def has_rings(mol):
    ring_functionals = ["[r5]", "[r6]", "[r7]", "[r8]"]
    ring_names = ["5-membered ring", "6-membered ring", "7-membered ring", "8-membered ring"]
    ring_groups = []
    ring_groups.append(ring_names)
    ring_groups.append(ring_functionals)
    ring_groups_in_mol = []
    for i in range(len(ring_groups[0])):
        pattern_smarts = Chem.MolFromSmarts(ring_groups[1][i])
        #print(type(mol))
        if mol.HasSubstructMatch(pattern_smarts) == True:
            ring_groups_in_mol.append(ring_groups[0][i])
    if len(ring_groups) != 0:
        return(ring_groups_in_mol)
    if len(ring_groups) == 0:
        no_list = ["N"]
        return(no_list)

def write_output(filepath, lines):
    f = open(filepath, "w")
    for line in lines:
        f.write(line)
        f.write('\n')
    f.close()   
    return()

def get_homo_lumo_from_xyz(xyz_filepath):
    atomic_numbers = {'C': 6, 'H': 1, 'O': 8, 'S': 16, 'F': 9, 'N': 7, 'Cl': 17}
    atoms = []

    with open(xyz_filepath, 'r') as xyz_file:
        for i, line in enumerate(xyz_file):
            # Skip the first two lines (header lines)
            if i < 2:
                continue
            
            # Split the line and extract the atomic symbol
            parts = line.split()
            if parts:  # Ensure the line isn't empty
                atom = parts[0].strip()
                # Verify that the atom is in the atomic_numbers dictionary
                if atom in atomic_numbers:
                    atoms.append(atom)
                else:
                    raise KeyError(f"Unknown atomic symbol '{atom}' in line: {line.strip()}")
    
    # Calculate the total atomic number and derive HOMO and LUMO indices
    atomic_num = sum(atomic_numbers[atom] for atom in atoms)
    homo_num = int(atomic_num / 2 - 1)
    lumo_num = homo_num + 1
    return homo_num, lumo_num

def SmilesToPDB(smiles_string, output_file):
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

def pdb_to_mol(pdb_filename):
    """
    Converts a PDB file to a molecule object.

    Parameters:
    - pdb_filename (str): Path to the PDB file.

    Returns:
    Chem.Mol or None: A molecule object if conversion is successful, else None.

    The function reads the content of the PDB file, converts it to a molecule object,
    and assigns chiral tags based on the molecular structure. If the conversion is
    unsuccessful or the PDB file is empty, None is returned.
    """
    with open(pdb_filename, 'r') as pdb_file:
        pdb_content = pdb_file.read()

    mol = Chem.MolFromPDBBlock(pdb_content)
    if mol is not None:
        AllChem.AssignAtomChiralTagsFromStructure(mol)

    return mol
