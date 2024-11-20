from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from modules.sw_basic_functions import vol_from_mol, vol_from_smiles

class orca_molecule:
    """
    A class to represent a molecule analyzed using ORCA computational chemistry software.

    Attributes:
    -----------
    name : str
        The name of the molecule.
    smiles : str
        The Simplified Molecular Input Line Entry System (SMILES) string representing the molecule's structure.
    mw : float
        The molecular weight of the molecule.
    peak_area : float
        The peak area from spectroscopy data, often related to the concentration or amount of the molecule.
    total_energy : float
        The total electronic energy of the molecule calculated by ORCA.
    homo_lumo_gap : float
        The energy gap between the highest occupied molecular orbital (HOMO) and the lowest unoccupied molecular orbital (LUMO).
    chemical_hardness : float
        A measure of the molecule's resistance to change in electron distribution; often calculated as half the HOMO-LUMO gap.
    dipole_moment : float
        The dipole moment of the molecule, representing the separation of charge within the molecule.
    polarizability : float
        The ability of the molecule to be polarized by an external electric field.
    volume : float
        The molecular volume, which can be related to the space the molecule occupies.

    Methods:
    --------
    __init__(self, name, smiles, mw, peak_area, total_energy, homo_lumo_gap, chemical_hardness, dipole_moment, polarizability, volume):
        Initializes the orca_molecule with the specified properties.
    """
    def __init__(self, name, smiles, mw, peak_area, total_energy, homo_lumo_gap, chemical_hardness, dipole_moment, polarizability, volume): # This is where the properties of the class are specified   
       self.name = name
       self.smiles = smiles
       self.mw = mw
       self.peak_area = peak_area
       self.total_energy = total_energy
       self.homo_lumo_gap = homo_lumo_gap
       self.chemical_hardness = chemical_hardness
       self.dipole_moment = dipole_moment
       self.polarizability = polarizability
       self.volume = volume

def csv_to_orca_class(csv_file):
    import csv
    with open(csv_file, 'r') as file:
      reader = csv.reader(file)
      molecules = []
      for molecule in reader:
          molecules.append(molecule) 
          
    # self, name, smiles, mw, core, r_groups, tracked_cuts, peak_area, total_energy, homo_lumo_gap, chemical_hardness, dipole_moment):  
      orca_molecules = []
      for i in range(len(molecules)):
          if i == 0:
              continue # skip header row
          else:
              molecule_name = molecules[i][0].strip() # Strip removes white space from the string
              molecular_weight = round(float(molecules[i][1]), 2)
              peak_area = molecules[i][2]
              smiles = molecules[i][3]    
              tot_energy = molecules[i][4]
              homo = molecules[i][5]
              lumo = molecules[i][6]
              homo_lumo_gap = molecules[i][7]
              chemical_hardness = molecules[i][8]
              dipole_moment = molecules[i][9]
              polarizability = molecules[i][10]
              # pulls the info from the csv for each mol
 
              mol = Chem.MolFromSmiles(smiles)            
              volume = vol_from_smiles(smiles)

               #Creating orca molecule object
              m1 = orca_molecule(molecule_name, smiles, molecular_weight, peak_area, tot_energy, homo_lumo_gap, chemical_hardness, dipole_moment, polarizability, volume)
               
              orca_molecules.append(m1)
    return(orca_molecules)



















