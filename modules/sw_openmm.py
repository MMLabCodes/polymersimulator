# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:50:12 2024

@author: danie
"""
import parmed as pm
import numpy as np
import pandas as pd
import csv
import itertools
import matplotlib.pyplot as plt
import shutil
import datetime
import time

from openmm import *
from openmm.app import *
from openmm.unit import *
from openmm.app.topology import *
#from openmmml import MLPotential

from openbabel import openbabel, pybel

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

def SmilesToPDB_GenerateRescode(smiles, name, directory, residue_code_csv):
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
    residue_codes = load_residue_codes(residue_code_csv)
    # Check if the name or smiles is already in the database
    existing_entry = find_existing_entry(residue_codes, name, smiles)
    if existing_entry:
        # Use existing residue code
        residue_code = existing_entry[2]  # Assuming code is the third column
    else:
        # Generate a unique 3-letter residue code excluding forbidden codes
        residue_code = generate_unique_residue_code(residue_codes, forbidden_codes)
        # Update the CSV file with the new entry
        update_residue_codes_csv(name, smiles, residue_code, residue_code_csv)
    # Replace default "UNL" codes in the PDB content - NOTE, "UNL" is also a forbidden code as it is the default code.
    pdb_filepath = directory + "/" + name + ".pdb"
    SmilesToPDB(smiles, pdb_filepath)    
    with open(pdb_filepath, "r") as pdb_file:
        lines = pdb_file.readlines()     
        for i, line in enumerate(lines):
            lines[i] = line.replace(" UNL", f" {residue_code}")   
    with open(pdb_filepath, "w") as pdb_file:              
        pdb_file.writelines(lines)
    return()

def load_residue_codes(residue_code_csv):
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

def generate_unique_residue_code(residue_codes, forbidden_codes=None):
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

def update_residue_codes_csv(name, smiles, residue_code, residue_code_csv):
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
    residue_codes = load_residue_codes(residue_code_csv)
    existing_entry = find_existing_entry(residue_codes, name, smiles)

    if not existing_entry:
        with open(residue_code_csv, "a", newline="") as csv_file:
            writer = csv.writer(csv_file)

            # Add a new entry to the CSV file
            writer.writerow([name, smiles, residue_code])

def find_existing_entry(residue_codes, name, smiles):
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
    for entry in residue_codes:
        if len(entry) >= 2 and (entry[0] == name or entry[1] == smiles):
            return entry
    return None

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

def get_rmsd(pdb_ref, pdb_prb):
    """
    Calculate the Root Mean Square Deviation (RMSD) between two molecular structures.

    Parameters
    ----------
    pdb_ref : str
        The path to the PDB file of the reference molecule.
    pdb_prb : str
        The path to the PDB file of the probe molecule.

    Returns
    -------
    float
        The RMSD between the reference and probe molecules after alignment.

    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    # Load your molecules
    prbMol = Chem.MolFromPDBFile(pdb_prb)
    refMol = Chem.MolFromPDBFile(pdb_ref)

    # Align the probe molecule to the reference molecule
    rmsd = rdMolAlign.AlignMol(prbMol, refMol)

    return rmsd



class TopologyML(Topology):
    """TopologyML holds the topology information needed by ANI force filed.
    There are errors with OpenMM-ML importing pdb files built using packmol due to the use of
    non standar residues.

    Parmed is used to load the pdb.
    """

    def __init__(self, input_pdb):
        """
        Parameters
        ----------
        self : TopologyML
            The new object created
        input_pdb : string
            Pdb file to load with parmed.
        """
        # Reading pdb with parmed
        self.pdb = pm.read_pdb(input_pdb)
        # Creating a PSF file with parmed
        self.pdb.write_psf(f'{input_pdb}.psf')
        # Createing OpenMM topology object
        self.new_topo = Topology()
        # Transforming to pandas dataframe the pdb
        self.pdb_pd = self.pdb.to_dataframe()

    def get_pm_chains(self):
        """Gets PDB chains from Parmed

        Return
        ------
        _chains : list
        """
        self._chains = self.pdb_pd['chain'].unique()  
        return self._chains
    
    def get_chain_obj(self, id):
        """Return the OpenMM chain object associated to the PDB chain readed.

        Parameters
        ----------
        id : string
            Chain identificator string.
        """
        for _c in self.new_topo.chains():
            if _c.id == id:
                return _c

    def create_chains(self, chains):
        """Creates OpenMM chains

        Parameters
        ----------
        chains : list
            Chains to be created.
        """
        for _c in chains:
            self.new_topo.addChain(_c)

    def get_pm_residues(self):
        """Gets residue list using parmed

        Returns
        -------
        _residues : list
            Residues list
        """
        self._residues = list(self.pdb_pd[['resname', 'chain','resnum']].groupby(['resname', 'chain', 'resnum']).groups)
        return self._residues

    def get_res_obj(self, chain, res_id):
        """ Returns the OpenMM residue object associated to the PDB residue
        created before.

        Parameters
        ----------
        chain : string
            Chain name in PDB file
        res_id : string

        Returns
        -------
        self._r : OpenMM residue object
        """
        self._c = self.get_chain_obj(chain)
        for self._r in self._c.residues():
            if self._r.id == int(res_id):
                return self._r

        
    def create_residue(self, resname, chain_name, resnum):
        """Create residue on Topology object given a set of properties.
        It is necessary to use the method `get_chain_name` to transform
        the chain_name to its associatd chain object created previously. 

        Parameters
        ----------
        resname : string
            Residue name given by parmed.
        resnum : string
            Residue number given by parmed.
        chain_name : string
            Chain name given by parmed

        Returns
        -------
        _residues : list
            Residues list
        """
        _chain_obj = self.get_chain_obj(chain_name)
        self.new_topo.addResidue(resname, _chain_obj, resnum)

    def create_residues(self, residues):
        for _c in self.new_topo.chains():
            for _r in residues:
                _resname = _r[0]
                _chain_name = _r[1]
                _resnum = _r[2]
                if _c.id == _chain_name:
                    self.create_residue(_resname, _chain_name, _resnum)

    def get_pm_atoms(self):
        """Gets lists of atoms from parmed.

        Returns
        -------
        atoms : list
            Atoms list of lists [[number, name, chain, resnum], ...] 
        """
        #self._atoms = self.pdb.atoms
        #return self._atoms
        return self.pdb
        
    def get_pm_num_atoms(self):
        self._n = 1
        for self._x in self.pdb:
            i = i + 1
        return i
    
    def create_topo(self):
        '''This method creates an OpenMM Topology object. The pdb file is parsed
        and translated to an OpenMM object.
        '''
        # Reading PDB chains with parmed.
        self._chains = self.get_pm_chains()
        # Reading PDB residues with parmed.
        self._residues = self.get_pm_residues()
        # Reading PDB atoms with parmed.
        self._atoms = self.get_pm_atoms()
        # Let's create
            # Chains
        self.create_chains(self._chains)
        self._n_chains = self.new_topo.getNumChains()
        print(f'Created {self._n_chains} chains')
            # Residues
        self.create_residues(self._residues)
        self._n_residues = self.new_topo.getNumResidues()
        print(f'Created {self._n_residues} residues')
        # Atoms
        # The chains and residues have been created in topology. Let's try with atoms.
        for self._atom in self._atoms:
            self.name = f'{self._atom.name}'
            self.element = Element.getBySymbol(self.name[0])
            res = self._atom.residue.name
            res_number = self._atom.residue.number
            number = self._atom.number
            chain = self._atom.residue.chain
            res_obj = self.get_res_obj(chain, res_number)
            self.new_topo.addAtom(self.name, self.element, res_obj)
        # Returning topology object
        return self.new_topo
    
class PositionsML():
    def __init__(self, input_pdb):
        """Extract the position from the pdb file. OpenMMM works in nanometer
        units, so, every coordinate system is converted to nanometer.

        Parameters
        ----------
        input_pdb : string
            Pdb file to load with parmed.
        
        Returns array
        -------
        coords: numpyt array with coordinates

        """
        # Reading pdb with parmed
        self.pdb = pm.read_pdb(input_pdb)
        # Createing OpenMM positions object
        self.positions = self.pdb.positions
        self.xyz = np.array(self.positions/nanometer)
        # The coordinates are translated to the origin.
        self.xyz[:,0] -= np.amin(self.xyz[:,0])
        self.xyz[:,1] -= np.amin(self.xyz[:,1])
        self.xyz[:,2] -= np.amin(self.xyz[:,2])
        self.coords = self.xyz * nanometer
    
    def get_pbc(self):
        """Generate the PBC vectors.
        
        """
        self.coords = self.coords/nanometer
        self.min_coord = np.min(self.coords, axis=0)
        self.max_coord = np.max(self.coords, axis=0)
        self.pbc_coord = self.max_coord - self.min_coord
        self.a = [self.pbc_coord[0], 0, 0]
        self.b = [0, self.pbc_coord[1], 0]
        self.c = [0, 0, self.pbc_coord[2]]
        return [self.a, self.b , self.c]

class DcdWriter():
    def __init__(self, prefix, freq):
        self.dcdReporter = app.DCDReporter(f'{prefix}.dcd', freq, enforcePeriodicBox = True)

class DataWriter():
    def __init__(self, prefix, freq, steps):
        self.stateDataReporter = app.StateDataReporter(f'{prefix}.txt', freq, totalSteps=steps,
step=True, time=True, speed=True, progress=True, elapsedTime=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, volume=True, density=True, separator=',')


class BuildSimulation():
    """
    Class for building and managing molecular simulations.

    Attributes:
        pressure (float): The default pressure for simulations, in atmospheres.
        temp (float): The default temperature for simulations, in Kelvin.
        timestep (float): The default timestep for simulations, in femtoseconds.
        friction_coeff (float): The default friction coefficient for Langevin dynamics simulations, in units of 1/picoseconds.
        total_steps (int): The default total number of steps for simulations.
        reporter_freq (int): The default frequency for reporting data during simulations.
        anneal_parameters (list): Parameters for simulated annealing, including start temp, max temp, cycles, holding steps, 
            and steps at each temperature.

    Methods:
        type_of_simulation(self): Determine the type of simulation.
        minimize_energy(self): Perform energy minimization of the system.
        minimize_energy_help(cls): Display help information for the minimize_energy method.
        anneal(self, simulation, start_temp=None, max_temp=None, cycles=None, holding_steps=None, steps_at_temp=None): 
            Perform simulated annealing on the provided simulation system.
        anneal_help(cls): Display help information for the anneal method.
        equilibrate(self, simulation, total_steps=None, temp=None, pressure=None): Equilibrate the provided simulation.
        equilibrate_help(cls): Display help information for the equilibrate method.
        production_run(self, simulation, total_steps=None, temp=None): Perform a production run simulation.
        production_run_help(cls): Display help information for the production_run method.
        __repr__(self): Representation of the simulation parameters.
        __str__(self): String representation of the simulation object.
        set_temperature(cls, temp): Set the temperature for simulations.
        set_pressure(cls, pressure): Set the pressure for simulations.
        set_timestep(cls, timestep): Set the timestep for simulations.
        set_friction_coeff(cls, friction_coeff): Set the friction coefficient for simulations.
        set_total_steps(cls, total_steps): Set the total number of steps for simulations.
        set_reporter_freq(cls, reporter_freq): Set the reporter frequency for simulations.
        set_anneal_parameters(cls, new_anneal_parameters): Set annealing parameters.
        graph_state_data(data_file): Graph the state data from a data file.

    Notes:
        This class provides a comprehensive interface for setting up, running, and managing molecular simulations.
        It includes methods for energy minimization, simulated annealing, equilibration, production runs, and visualization of data.

    Child Classes (this class is used through these child classes, do not try to call BuildSimulaiton on its own unless you plan to only utilise static methods):
        - ANISimulation: A class representing ANI (Atomic Neural Network Interaction) molecular dynamics simulations.
        - AmberSimulation: A class representing AMBER molecular dynamics simulations.
    """ 
    pressure = 1
    temp = 300
    timestep = 2.0
    friction_coeff = 1.0
    total_steps = 1000
    reporter_freq = 1000
    nonbondedcutoff = 1.0
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    # Start temp, max temp, cycles, holding_steps, steps at each temp
    anneal_parameters = [300, 700, 5, 3000, 100]
    minimized_only = None # This will change to True/False and will determine how periodic box vectors are set
    
    def __init__(self, manager, filename):
        """
        Initializes a BuildSimulation object with specified directories and filename 

        Args:
            directories (str): The directories where simulation files are stored.
            filename (str): The name of the simulation file.

        Attributes:
            filename (str): The name of the simulation file.
            directories (str): The directories where simulation files are stored.
            output_dir (str): The directory for simulation output.
            log_info (dict): Information log for different stages of the simulation.
            log_csv (str): The path to the CSV log file.

        Notes:
            If the output directory does not exist, it will be created.
        """
        self.manager = manager
        self.filename = filename
        self.output_dir = os.path.join(self.manager.systems_dir, self.filename, self.timestamp)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.log_info = {'Minimization': {'Temperature':None,'Time taken': None},
                            'Annealing': {'Time taken': None, 'Simulation time': None, 'Start temp': None, 'Target temp': None, 'Cycles': None, 'Steps at plateaus': None, 'Steps at incremental temps': None, 'Timestep': None},
                            'Equilibration': {'Time taken': None, 'Simulation time': None, 'Temperature': None, 'Pressure': None, 'Timestep': None},
                            'Production': {'Time taken': None, 'Simulation time': None, 'Temperature': None, 'Timestep': None},
                        'Thermal Ramp' : {'Time taken':None, 'Simulation time':None, 'Start temp':None, 'Target temp':None, 'Quench Rate':None, 'Steps at incremental temps':None, 'Timestep':None, 'Ensemble':None, 'Method':None}}
        self.log_csv = os.path.join(self.output_dir, (self.filename + "_" + self.timestamp + "_log.csv"))
                             
    def type_of_simulation(self):
        """
        Determine the type of simulation.
        
        Returns:
            str: A string indicating the type of simulation, either "AMB" for AmberSimulation or "ANI" for ANISimulation.

        Notes:
            This method inspects the type of the current object to determine the type of simulation.
            If the object is an instance of AmberSimulation, it returns "AMB".
            If the object is an instance of ANISimulation, it returns "ANI".
            If the object is neither, it returns a message prompting the user to specify the type of simulation by calling 
            either AmberSimulation or ANISimulation classes.
            
            This function should never need to be called as we can we can the type of our simulation object with: repr(simulation_object)
        """
        if isinstance(self, AmberSimulation):
            return "AMB"
        if isinstance(self, ANISimulation):
            return "ANI"
        else:
            return "Please specify type of simulation by calling <AmberSimulation> OR <ANISimulation> classes"
    
    def minimize_energy(self):
        """
        Function to perform energy minimization of the system using Langevin dynamics.
        
        USAGE:
            minimized_sim = simulation_object.minimize_energy()

        Returns:
            minimized_simulation_object
        """
        min_start_time = time.time()
        
        integrator = LangevinIntegrator(self.temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator)
            simulation.context.setPositions(self.amb_coordinates.positions)
            
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*unit.nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)
            simulation.context.setPositions(self.ani_coordinates)
        
        simulation.minimizeEnergy()

        state = simulation.context.getState(getPositions=True, getEnergy=True) 
        
        self.min_pdbname = os.path.join(self.output_dir, ("min_" + self.filename + ".pdb"))
        with open(self.min_pdbname, 'w') as output:
            PDBFile.writeFile(simulation.topology, state.getPositions(), output)
        
        min_end_time = time.time()  
        time_taken = min_end_time - min_start_time
        self.log_info['Minimization']['Time taken'] = time_taken
        self.log_info['Minimization']['Temperature'] = self.temp
        return(simulation)
     
    @classmethod
    def minimize_energy_help(cls):
        """Display help information for the minimize_energy method."""
        print(cls.minimize_energy.__doc__)    
     
    def anneal(self, simulation, start_temp=None, max_temp=None, cycles=None, holding_steps=None, steps_at_temp=None):
        """
        Function to perform simulated annealing on the provided simulation system.
        
        USAGE: 
            annealed_sim = sim.anneal(simulation, start_temp, max_temp, cycles, holding_steps, steps_at_temp)            
        
        Recommended USAGE:
            annealed_sim = sim.anneal(simulation)
            
            where annealing parameters are set with the following function:
                
            simulation_object.set_anneal_parameters([start_temp, max_temp, cycles, holding_steps, steps_at_each_temp])
        
        Args:
            simulation (app.Simulation): The simulation object to perform annealing on.
            
            start_temp (float, optional): The starting temperature for annealing in Kelvin. Defaults to None, 
                in which case the value is fetched from self.anneal_parameters[0].
                
            max_temp (float, optional): The maximum temperature for annealing in Kelvin. Defaults to None, 
                in which case the value is fetched from self.anneal_parameters[1].
                
            cycles (int, optional): The number of annealing cycles to perform. Defaults to None, 
                in which case the value is fetched from self.anneal_parameters[2].
                
            holding_steps (int, optional): The number of steps to hold the system at each temperature. 
                Defaults to None, in which case the value is fetched from self.anneal_parameters[3].
                
            steps_at_temp (int, optional): The number of steps to perform at each temperature. 
                Defaults to None, in which case the value is fetched from self.anneal_parameters[4].

        Returns:
            A tuple containing the simulation object after annealing and the filename of the data file generated.
                in the following format:
        
            (simulation_object, path_to_data_file)
                    
        Notes:
            This method performs simulated annealing on the provided simulation object. It initializes the system with the 
        provided initial conditions and runs annealing cycles adjusting the temperature according to the specified parameters.
        The system's state is updated throughout the annealing process, and reporters are set up to record trajectory and 
        simulation data.
        """
        anneal_start_time = time.time()
        
        if start_temp is None:
            start_temp = self.anneal_parameters[0]
        if max_temp is None:
            max_temp = self.anneal_parameters[1]
        if cycles is None:
            cycles = self.anneal_parameters[2]
        if holding_steps is None:
            holding_steps = self.anneal_parameters[3]
        if steps_at_temp is None:
            steps_at_temp = self.anneal_parameters[4]
         
        # Extract positional info
        state = simulation.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors() # Obtain periodc box vectors of the previous step
        
        # Set up integrator
        integrator = LangevinIntegrator(start_temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*nanometers, constraints=app.HBonds)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)

        # Update the xyz of each atom
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        simulation.context.setPositions(xyz)
        
        # Total up all steps of the equilibration so the rpeorters record properly
        total_steps = (((max_temp-start_temp)*1000 + holding_steps)*2)*cycles
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for        
        #output_pdbname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_anneal.pdb"))
        output_pdbname = os.path.join(self.output_dir, (self.filename + "_anneal_" + str(self.timestamp) + ".pdb" ))
        simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
        
        # DCD trajectory
        #output_dcdname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_anneal"))
        output_dcdname = os.path.join(self.output_dir, (self.filename + "_anneal_traj_" + str(self.timestamp)))
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        #output_dataname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_anneal_data"))
        output_dataname = os.path.join(self.output_dir, (self.filename + "_anneal_data_" + str(self.timestamp)))
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter)
        
        increments = int(max_temp-start_temp)
        
        for i in range(cycles):
            for j in range(increments):
                integrator.setTemperature(start_temp+j)
                simulation.step(steps_at_temp)
            
            integrator.setTemperature(max_temp)
            simulation.step(holding_steps)
            
            for j in range(increments):
                integrator.setTemperature(max_temp-j)
                simulation.step(steps_at_temp)
                
            integrator.setTemperature(start_temp)
            simulation.step(holding_steps)

        # Now copy all the files into the system directory
        #destination_dir = os.path.join(directories.systems_dir, self.filename)
        #shutil.move(output_pdbname, destination_dir)
        #output_dcdname = self.filename +  "_anneal_traj_" + self.timestamp + ".dcd"
        #shutil.move(output_dcdname, destination_dir)
        #output_dataname = self.filename +  "_anneal_data_" + self.timestamp + ".txt"
        #shutil.move(output_dataname, destination_dir)
        #return(simulation, os.path.join(destination_dir, output_dataname))

        anneal_end_time = time.time()  
        time_taken = anneal_end_time - anneal_start_time

        self.log_info['Annealing']['Time taken'] = time_taken
        self.log_info['Annealing']['Simulation time'] = total_steps*self.timestep
        self.log_info['Annealing']['Start temp'] = start_temp
        self.log_info['Annealing']['Target temp'] = max_temp
        self.log_info['Annealing']['Cycles'] = cycles
        self.log_info['Annealing']['Steps at plateaus'] = holding_steps
        self.log_info['Annealing']['Steps at incremental temps'] = steps_at_temp
        self.log_info['Annealing']['Timestep'] = self.timestep
        return(simulation, (output_dataname + ".txt"))
    
    @classmethod
    def anneal_help(cls):
        """Display help information for the anneal method."""
        print(cls.anneal.__doc__)    
    
    def equilibrate(self, simulation, total_steps=None, temp=None, pressure=None):
        """
        Function to equilibrate the provided simulation to reach a specified temperature and pressure.
        
        USAGE: 
            equilibrated_sim = sim.anneal(simulation, total_steps, temp, pressure)            
        
        Recommended USAGE:
            equilibrated_sim = sim.anneal(simulation)
            
            where annealing parameters are set with the following functions:
                
            simulation_object.set_temperature(temperature)
            simulation_object.set_pressure(pressure)
            simulation_object.set_total_steps(total_steps)
            
        Args:
            simulation (app.Simulation): The simulation object to equilibrate.
            
            total_steps (int, optional): The total number of steps to run for equilibration. Defaults to None, 
                in which case the value is fetched from self.total_steps.
                
            temp (float, optional): The target temperature for equilibration in Kelvin. Defaults to None, 
                in which case the value is fetched from self.temp.
                
            pressure (float, optional): The target pressure for equilibration in atmospheres. Defaults to None, 
                in which case the value is fetched from self.pressure.

        Returns:
            A tuple containing the simulation object after equilibration and the filename of the data file generated
                in the following format:
                    
            (simulation_object, path_to_data_file)

        Notes:
            This method performs equilibration on the provided simulation object to reach the specified temperature and pressure.
            It initializes the system with the provided initial conditions and runs the simulation for the specified number 
            of steps. The system's state is updated throughout the equilibration process, and reporters are set up to record 
            trajectory and simulation data.
        """
        equili_start_time = time.time()
        
        if total_steps is None:
            total_steps = self.total_steps
        if temp is None:
            temp = self.temp
        if pressure is None:
            pressure = self.pressure
            
        # Extract positional info
        state = simulation.context.getState(getPositions=True, getEnergy=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors()
        
        # Set up integrator and barostat
        barostat = MonteCarloBarostat((pressure*atmosphere), (temp*kelvin)) # Define barostat (pressure, temp)
        integrator = LangevinIntegrator(temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*nanometers, constraints=app.HBonds)
            system.addForce(barostat)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*unit.nanometers, constraints=app.HBonds)
            system.addForce(barostat)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)

        # We set the box vectors with the output from the the previous simulation
        #if self.minimized_only == False:       
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        
        # Update positional info
        simulation.context.setPositions(xyz)
             
        # Set initial velocities
        simulation.context.setVelocitiesToTemperature(temp*kelvin)
        
        # Set up reporters
        #output_pdbname = os.path.join(directories.systems_dir, self.filename, (self.filename +  "_" + str(pressure) + "_atm_traj.pdb"))
        output_pdbname = os.path.join(self.output_dir, (self.filename +  "_" + str(pressure) + "_atm_traj_" + str(self.timestamp) + ".pdb"))
        simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
    
        # DCD trajectory
        #output_dcdname = os.path.join(directories.systems_dir, self.filename, (self.filename +  "_" + str(pressure) + "_atm_traj.dcd"))
        output_dcdname = os.path.join(self.output_dir, (self.filename +  "_" + str(pressure) + "_atm_traj_" + str(self.timestamp)))
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        #output_dataname = os.path.join(directories.systems_dir, self.filename, (self.filename +  "_" + str(pressure) + "_atm_data"))
        output_dataname = os.path.join(self.output_dir, (self.filename +  "_" + str(pressure) + "_atm_data_" + str(self.timestamp)))
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter)
        simulation.step(total_steps)       

        # Now copy all the files into the system directory
        #destination_dir = os.path.join(directories.systems_dir, self.filename)
        #shutil.move(output_pdbname, destination_dir)
        #output_dcdname = self.filename +  "_" + str(pressure) + "_atm_traj_" + str(self.timestamp) + ".dcd"
        #shutil.move(output_dcdname, destination_dir)
        #output_dataname = self.filename +  "_" + str(pressure) + "_atm_data_" + str(self.timestamp) + ".txt"
        #shutil.move(output_dataname, destination_dir)
        equili_end_time = time.time()  
        time_taken = equili_end_time - equili_start_time

        self.log_info['Equilibration']['Time taken'] = time_taken
        self.log_info['Equilibration']['Simulation time'] = total_steps * self.timestep
        self.log_info['Equilibration']['Temperature'] = temp
        self.log_info['Equilibration']['Pressure'] = pressure
        self.log_info['Equilibration']['Timestep'] = self.timestep
        return(simulation, (output_dataname + ".txt"))
    
    @classmethod
    def equilibrate_help(cls):
        """Display help information for the equilibrate method."""
        print(cls.equilibrate.__doc__) 
        
    def production_run(self, simulation, total_steps=None, temp=None):
        """
        Function to perform a production run simulation with the provided parameters.
        
        USAGE: 
            production_run_sim = sim.anneal(simulation, total_steps, temp)            
        
        Recommended USAGE:
            production_run_sim = sim.anneal(simulation)
            
            where annealing parameters are set with the following functions:
                
            simulation_object.set_temperature(temperature)
            simulation_object.set_total_steps(total_steps)
            
        Args:
            simulation (app.Simulation): The simulation object to run the production simulation on.
            
            total_steps (int, optional): The total number of steps to run for the production simulation. 
                Defaults to None, in which case the value is fetched from self.total_steps.
            
            temp (float, optional): The temperature for the production simulation in Kelvin. 
                Defaults to None, in which case the value is fetched from self.temp.

        Returns:
            A tuple containing the simulation object after the production run and the filename of the data file generated.
                in the following format:
                
            (simulation_object, path_to_data_file)

        Notes:
            This method performs a production run simulation with the provided simulation object and parameters.
            It initializes the system with the provided initial conditions and runs the simulation for the specified 
            number of steps. The system's state is updated throughout the production run, and reporters are set up to 
            record trajectory and simulation data.
        """
        prod_start_time = time.time()
        if total_steps is None:
            total_steps = self.total_steps
        if temp is None:
            temp = self.temp
            
        # Extract positional info
        state = simulation.context.getState(getPositions=True, getEnergy=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors() # Obtain periodc box vectors of the previous step
        
        # Set up integrator
        integrator = LangevinIntegrator(temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*nanometers, constraints=app.HBonds)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)  
        
        # Update positional info
        simulation.context.setPositions(xyz)
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        #   visualisation of the system and coloring specific residues 
        #output_pdbname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_prod_traj.pdb"))
        output_pdbname = os.path.join(self.output_dir, (self.filename + "_prod_traj_" + str(self.timestamp) + ".pdb"))
        simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
        
        # DCD trajectory
        #output_dcdname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_prod_traj"))
        output_dcdname = os.path.join(self.output_dir, (self.filename + "_prod_traj_" + str(self.timestamp)))
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        #output_dataname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_prod_data"))
        output_dataname = os.path.join(self.output_dir, (self.filename + "_prod_data_" + str(self.timestamp)))
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter) 
        simulation.step(total_steps)

        # Now copy all the files into the system directory
        #destination_dir = os.path.join(directories.systems_dir, self.filename)
        #shutil.move(output_pdbname, destination_dir)
        #output_dcdname = self.filename +  "_prod_traj_" + str(self.timestamp) + ".dcd"
        #shutil.move(output_dcdname, destination_dir)
        #output_dataname = self.filename +  "_prod_data_" + str(self.timestamp) + ".txt"
        #shutil.move(output_dataname, destination_dir)
        prod_end_time = time.time()
        time_taken = prod_end_time - prod_start_time
        self.log_info['Production']['Time taken'] = time_taken
        self.log_info['Production']['Simulation time'] = total_steps * self.timestep
        self.log_info['Production']['Temperature'] = temp
        self.log_info['Production']['Timestep'] = self.timestep
        return(simulation, (output_dataname + ".txt"))
        
    def thermal_ramp(self, simulation, heating, quench_rate, ensemble, start_temp=None, max_temp=None, total_steps=None, pressure=None):
        # Heating is a boolean - true for heating, false for cooling
        # quench rate is an integer. i.e. 10 or 20 K (or smth else) and is the temperature steps the system is heated/cooled in
        # Ensemble is "NVT" or "NPT"
        if ensemble == "NVT" or ensemble == "NPT":
            pass
        else:
            print("Please specify 'NVT' or 'NPT' ensemble for the thermal ramp")
            return()
            
        thermal_ramp_start_time = time.time()

        if start_temp is None:
            start_temp = self.anneal_parameters[0]
        if max_temp is None:
            max_temp = self.anneal_parameters[1]
        if pressure is None:
             pressure = self.pressure
        if total_steps is None:
            total_steps = self.total_steps

        # Extract positional info
        state = simulation.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors() # Obtain periodc box vectors of the previous step
        
        # Set up integrator
        integrator = LangevinIntegrator(start_temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*nanometers, constraints=app.HBonds)
            if ensemble == "NPT":
                barostat = MonteCarloBarostat((pressure*atmosphere), ((start_temp if heating==True else max_temp)*kelvin))
                system.addForce(barostat)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff1.0*nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)

        # Update the xyz of each atom
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        simulation.context.setPositions(xyz)

        method = "heat" if heating==True else "cool"
        output_filename = self.filename + "_temp_ramp_" + method + f"_{str(start_temp)}_{str(max_temp)}_" + str(self.timestamp)
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        output_pdbname = os.path.join(self.output_dir, (output_filename + ".pdb" ))
        simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
        
        # DCD trajectory
        output_dcdname = os.path.join(self.output_dir, output_filename)
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        output_dataname = os.path.join(self.output_dir, output_filename)
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter)
        
        incremental_temps = np.arange(start_temp, max_temp + quench_rate, quench_rate).tolist()
     
        if heating is True:
            pass # List will be low to high
        if heating is False:
            incremental_temps.reverse() # List now high to low for cooling
         
        steps_at_increment = int(total_steps/len(incremental_temps))

        for i in range(len(incremental_temps)):
            integrator.setTemperature(incremental_temps[i])
            simulation.step(steps_at_increment)

        thermal_ramp_end_time = time.time()
        time_taken = thermal_ramp_start_time - thermal_ramp_end_time

        self.log_info['Thermal Ramp']['Time taken'] = time_taken
        self.log_info['Thermal Ramp']['Simulation time'] = total_steps*self.timestep
        self.log_info['Thermal Ramp']['Start temp'] = start_temp
        self.log_info['Thermal Ramp']['Target temp'] = max_temp
        self.log_info['Thermal Ramp']['Quench Rate'] = quench_rate
        self.log_info['Thermal Ramp']['Steps at incremental temps'] = steps_at_increment
        self.log_info['Thermal Ramp']['Timestep'] = self.timestep
        self.log_info['Thermal Ramp']['Ensemble'] = ensemble
        self.log_info['Thermal Ramp']['Method'] = "heating" if heating==True else "cooling" 

        return(simulation, (output_dataname + ".txt"))


    
    @classmethod
    def production_run_help(cls):
        """Display help information for the production run method."""
        print(cls.production_run.__doc__)
        
    def __repr__(self):
        print("Simulation parameters: ('{}', '{}', '{}, {}, {}, {}')".format(self.pressure, self.temp, self.timestep, self.friction_coeff, self.total_steps, self.reporter_freq))
        return("Simulation parameters given in the following format: ('{}', '{}', '{}, {}, {}, {}')".format("pressure", "temperature", "timestep", "friction coefficient", "total steps", "reporter freqeuncy"))
    
    def __str__(self):
        return 'Simulation object of - {}'.format( self.filename)
    
    # All class methods here change a variable in the following way
    # AmberSimulation.set_temperature(400)
    # AmberSimulation.temp

    @classmethod
    def display_start_time(cls):
        print("Simulation initiated at: ", cls.timestamp)
    
    @classmethod
    def set_temperature(cls, temp): 
        """
        Class method to set the temperature for simulations.

        Args:
            cls: The class itself.
            temp (float): The temperature to set, in Kelvin.

        Returns:
            None

        Notes:
            This method sets the temperature class attribute to the specified value and prints a confirmation message.
        """
        cls.temp = temp
        print("Temperature set to: ", str(temp), "kelvin")
             
    @classmethod
    def set_pressure(cls, pressure):
        """
        Class method to set the pressure for simulations.

        Args:
            cls: The class itself.
            pressure (float): The pressure to set, in Atmospheres.

        Returns:
            None

        Notes:
            This method sets the pressure class attribute to the specified value and prints a confirmation message.
        """
        cls.raise_pressure = pressure
        print("Pressure set to: ", str(pressure), " atmospheres")
        
    @classmethod
    def set_timestep(cls, timestep): 
        """
        Class method to set the timestep for simulations.

        Args:
            cls: The class itself.
            timestep (float): The timestep to set, in Femtoseconds.

        Returns:
            None

        Notes:
            This method sets the timestep class attribute to the specified value and prints a confirmation message.
        """
        cls.timestep = timestep
        print("Timestep set to: ", str(timestep))
        
    @classmethod
    def set_friction_coeff(cls, friction_coeff):
        """
        Class method to set the friction_coeff for simulations.

        Args:
            cls: The class itself.
            friction_coeff (float): The friction_coeff to set, in Picoseconds.

        Returns:
            None

        Notes:
            This method sets the friction_coeff class attribute to the specified value and prints a confirmation message.
        """
        cls.friction_coeff = friction_coeff
        print("Friction coeffiecent set to: ", str(friction_coeff))
    
    @classmethod
    def set_total_steps(cls, total_steps):
        """
        Class method to set the total_steps for simulations.

        Args:
            cls: The class itself.
            total_steps (int): The total_steps to set, as an integer.

        Returns:
            None

        Notes:
            This method sets the total_steps class attribute to the specified value and prints a confirmation message.
        """
        cls.total_steps = total_steps
        print("Total steps for simulation set to: ", str(total_steps))
        
    @classmethod
    def set_reporter_freq(cls, reporter_freq):
        """
        Class method to set the reporter_freq for simulations.

        Args:
            cls: The class itself.
            reporter_freq (int): The reporter_freq to set, as an integer.

        Returns:
            None

        Notes:
            This method sets the reporter_freq class attribute to the specified value and prints a confirmation message.
        """
        cls.reporter_freq = reporter_freq
        print("Reporter frequency set to every: ", reporter_freq, " steps")

    @classmethod
    def set_nonbondedcutoff(cls, nonbondedcutoff):
        if nonbondedcutoff is not type(float):
            print("Please pass a float to this method to set a new nonbondedcutoff")
            print("Example: simulation.set_nonbondedcutoff(5.0)")
            return(none)
        cls.nonbondedcutoff = nonbondedcutoff
        
    @classmethod
    def set_anneal_parameters(cls, new_anneal_parameters): 
        """
        Class method to set annealing parameters.

        Args:
            cls: The class itself.
            new_anneal_parameters (list): List of annealing parameters in the format [start_temp, max_temp, cycles, holding_steps, steps_at_temp].

        Returns:
            None

        Raises:
            ValueError: If the length of new_anneal_parameters does not match the expected length.

        Notes:
            This method sets the annealing parameters class attribute to the specified list. 
            It prints a confirmation message with the provided parameters.
        """
        if len(new_anneal_parameters) != len(cls.anneal_parameters):
            format_str = "Expected format: [start_temp, max_temp, cycles, holding_steps, steps_at_temp]"
            raise ValueError(f"Invalid parameters provided. {format_str}")
        else: 
            # new_anneal_parameters = [start_temp, max_temp, cycles, holding_steps]
            cls.anneal_parameters = new_anneal_parameters
            print("Anneal parameters set.")
            print("Starting temperature is: ", str(new_anneal_parameters[0]))
            print("Target temperature is: ", str(new_anneal_parameters[1]))
            print("Number of annealing cycles is: ", str(new_anneal_parameters[2]))
            print("Steps at target/start temperature is: ", str(new_anneal_parameters[3]))
            print("Steps at each incremental temperature is: ", str(new_anneal_parameters[4]))

    @classmethod
    def set_anneal_parameters_help(cls):
        """Display help information for setting annealing parameters."""
        print(cls.set_anneal_parameters.__doc__)
           
    @staticmethod  # Dont take instance or class - can run this without an instance
    def graph_state_data(data_file):
        """
        Function to generate a plot from simulation data and save it as a PNG file.
        
        USAGE:
            The variable used for this plot is the second element of the tuple coming from the following method:
                **anneal**
                **equilibrate**
                **production_run**
                
            sim.graph_state_data(returned_simulation_variable)
            
            An example of this returned_simulation_varaible is:
                sim.anneal(minimized_sim)[0]
        
        Args:
            data_file (str): The path to the CSV file containing simulation data.
            
        Notes:
            This method reads simulation data from a CSV file, creates plots for each column of data
            against time, and saves the resulting plot as a PNG file. The CSV file is expected to have
            a "Time (ps)" column and numerical data columns representing different states of the
            simulation over time.

            The generated plot will contain multiple subplots, each showing a different state variable
            plotted against time. The number of rows and columns of subplots is determined automatically
            based on the number of data columns, with a maximum of 2 columns per row. If there are more
            data columns than can fit in the specified layout, some subplots may be omitted.

            The PNG file containing the plot will be saved in the same directory as the data file with
            the same name, but with the extension changed to ".png".
        """
        png_file_name = (data_file.split(".")[0]) + ".png"
        
        data_file = data_file
        df = pd.read_csv(data_file, delimiter=',')

        # Exclude the last two columns
        columns_to_plot = df.columns[3:-2]

        # Calculate the number of rows and columns needed for subplots
        num_rows = (len(columns_to_plot) + 1) // 2  # +1 to ensure an extra row if there's an odd number of columns
        num_cols = min(2, len(columns_to_plot))  # Maximum of 2 columns in each row

        # Set up subplots in the calculated arrangement
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 4 * num_rows))

        # Plot each column against "Time (ps)" in a separate subplot
        for i, column in enumerate(columns_to_plot):
            row = i // num_cols
            col = i % num_cols
            axes[row, col].plot(df["Time (ps)"], df[column])
            axes[row, col].set_title(column)
            axes[row, col].set_xlabel("Time (ps)")
            axes[row, col].set_ylabel(column)
            axes[row, col].grid(True)

        # Remove empty subplots if there are not enough columns to fill all slots
        for i in range(len(columns_to_plot), num_rows * num_cols):
            row = i // num_cols
            col = i % num_cols
            fig.delaxes(axes[row, col])

        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust the layout to prevent clipping of the suptitle
        plt.savefig(png_file_name, dpi=600)
        plt.show()
        
    @classmethod
    def graph_state_data_help(cls):
        """Display help information for the method to plot state data."""
        print(cls.graph_state_data.__doc__)

    def write_log_csv(self):
        """
        Write the log information to a CSV file.

        Returns:
            None

        Notes:
            This method collects all unique keys across all sections of log_info.
            It then writes the log_info to a CSV file with each section as a row and the keys as columns.
        """
        # Get all unique keys across all sections
        all_keys = set()
        for info in self.log_info.values():
            all_keys.update(info.keys())

        # Write log_info to a CSV file
        with open(self.log_csv, 'w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['Section'] + list(all_keys))
            writer.writeheader()
            for section, info in self.log_info.items():
                writer.writerow({'Section': section, **info})
        
class AmberSimulation(BuildSimulation):
    """
    A class representing an AMBER molecular dynamics simulation.

    Attributes:
        directories (str): A directories obejct generated with sw_directories which manages filepaths during simulations.
        topology_file (str): A string representing the path to the topology file.
        coordinates_file (str): A string representing the path to the coordinates file.

    Methods:
        __init__(directories, topology_file, coordinates_file): 
            Initializes an AmberSimulation object with specified directories, topology file, and coordinates file.
        __str__():
            Returns a string representation of the AmberSimulation object.
    """
    def __new__(cls, *args, **kwargs):
        # Check number of arguments provided #
        # The __new__ method essentially creates a dummy instance before actually creating the real instance # 
        if len(args) != 3:
            raise TypeError(
                f"AmberSimulation expected 3 arguments, but {len(args)} were given. "
                "Usage: 'sim = AmberSimulation(manager, topology_file, coordinates_file)'"
            )
        # Call __new__ to continue object creation if the argument count is correct
        return super().__new__(cls)   
        
    def __init__(self, manager, topology_file, coordinates_file): 
        self.manger = manager
        self.filename = os.path.basename(topology_file).split('.')[0]
        
        # Inherit attributes defined in the parent class and pass directories to the BuildSimulation constructor
        super().__init__(manager, self.filename)
        
        self.amb_coordinates = app.AmberInpcrdFile(coordinates_file)
        self.amb_topology = app.AmberPrmtopFile(topology_file, periodicBoxVectors=self.amb_coordinates.boxVectors)
        
        # Set pbc into topology
        #self.positionsML = PositionsML(pdb_file)
        #self.pbc = self.positionsML.get_pbc()
        #self.amb_topology.topology.setPeriodicBoxVectors(self.pbc)

    def __str__(self):
        return 'Amber simulation object of - {}'.format(self.filename)

class ANISimulation(BuildSimulation):
    """
    A class representing an ANI (Atomic Neural Network Interaction) molecular dynamics simulation.

    Attributes:
        potential (str): A string representing the ANI potential used in the simulation.
        directories (str): A string representing the directories where simulation files are stored.
        input_file (str): A string representing the path to the input file.

    Methods:
        __init__(directories, input_file): 
            Initializes an ANISimulation object with specified directories and input file, sets up ANI coordinates, potential, and topology.
        __str__():
            Returns a string representation of the ANISimulation object.
        set_potential(potential):
            Class method to set the ANI potential for all instances of the class.
    """
    potential = 'ani2x'
    
    def __init__(self, directories, input_file):
        self.filename = os.path.basename(input_file).split('.')[0]
        
        # Inherit attributes defined in the parent class and pass directories to the BuildSimulation constructor
        super().__init__(directories, self.filename)
        
        self.positionsML = PositionsML(input_file)
        self.ani_coordinates = self.positionsML.coords
        self.pbc = self.positionsML.get_pbc()
        self.potential = MLPotential(self.potential)
        self.ani_topology = ((TopologyML(input_file)).create_topo()).setPeriodicBoxVectors(self.pbc)

    
    def __str__(self):
        return 'ANI simulation object of - {}'.format(self.filename)
    
    @classmethod
    def set_potential(cls, potential):
        # No need for this yet, but can be used to swap to ANI 1 potential
        cls.potential = potential
        print("Potential set to: ", str(potential))


'''
['number', 'name', 'type', 'atomic_number', 'charge', 'mass', 'nb_idx', 'solvent_radius', 'screen', 'occupancy', 'bfactor', 'altloc', 'tree', 'join', 'irotat', 'rmin', 'epsilon', 'rmin_14', 'epsilon_14', 'resname', 'resid', 'resnum', 'chain', 'segid', 'xx', 'xy', 'xz']
'''
