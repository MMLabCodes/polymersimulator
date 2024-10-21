from rdkit import Chem
from rdkit.Chem import Draw, AllChem

def vol_from_smiles(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol) 
    return volume

def vol_from_mol(mol):
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol)   
    return volume

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