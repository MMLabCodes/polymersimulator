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
    atomic_numbers = {'C': 6, 'H': 1, 'O': 8, 'S':16, 'F':9, 'N':14, 'Cl':17}
    atoms = []

    with open(xyz_filepath, 'r') as xyz_file:
        for i, line in enumerate(xyz_file):
            if i == 0:
                continue
            if i == 1:
                continue
            else: 
                atom = (line.split(" ")[0]).strip()
                atoms.append(atom)
                
    atomic_num = sum(atomic_numbers[atom] for atom in atoms)
    homo_num = int(atomic_num/2 - 1)
    lumo_num = homo_num + 1
    return(homo_num, lumo_num)