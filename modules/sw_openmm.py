# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:50:12 2024

@author: danie
"""
import parmed as pm
from parmed import load_file
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
        """
        Calculate and return the number of atoms in the `pdb` attribute.

        The function iterates through the `pdb` attribute (assumed to be an iterable 
        containing atomic data) and counts the number of elements (atoms).

        Returns:
            int: The number of atoms in the `pdb` attribute.
        """
        count = 0
        for _ in self.pdb:
            count += 1
        return count

    
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

class DcdWriter:
    """
    A class to handle writing DCD (Dynamics Coordinate Data) files for molecular simulations.

    Attributes:
        dcdReporter (app.DCDReporter): An instance of DCDReporter to write the DCD file with specified parameters.

    Methods:
        __init__(prefix: str, freq: int):
            Initializes the DcdWriter with the given prefix and frequency for writing frames.
    """

    def __init__(self, prefix, freq):
        """
        Initialize a DcdWriter instance.

        Args:
            prefix (str): The prefix for the output DCD file name.
            freq (int): The frequency (in steps) at which frames are written to the DCD file.

        The output file will be named `{prefix}.dcd`. The DCDReporter enforces the periodic box
        constraints for the simulation.
        """
        self.dcdReporter = app.DCDReporter(f'{prefix}.dcd', freq, enforcePeriodicBox=True)


class DataWriter:
    """
    A class to handle writing simulation state data to a text file.

    Attributes:
        stateDataReporter (app.StateDataReporter): An instance of StateDataReporter to log simulation data.

    Methods:
        __init__(prefix: str, freq: int, steps: int):
            Initializes the DataWriter with the given prefix, frequency, and total number of steps.
    """

    def __init__(self, prefix, freq, steps):
        """
        Initialize a DataWriter instance.

        Args:
            prefix (str): The prefix for the output text file name.
            freq (int): The frequency (in steps) at which data is written to the text file.
            steps (int): The total number of simulation steps.

        The output file will be named `{prefix}.txt`. The StateDataReporter logs various
        simulation parameters, including step count, time, speed, progress, elapsed time, 
        total energy, kinetic energy, potential energy, temperature, volume, and density.
        The data is separated by commas.
        """
        self.stateDataReporter = app.StateDataReporter(
            f'{prefix}.txt', freq, totalSteps=steps,
            step=True, time=True, speed=True, progress=True, elapsedTime=True,
            totalEnergy=True, kineticEnergy=True, potentialEnergy=True,
            temperature=True, volume=True, density=True, separator=','
        )



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
    savepdb_traj = False
    pressure = 1
    temp = 300
    min_temp = 0
    timestep = 2.0
    friction_coeff = 1.0
    total_steps = 1000
    reporter_freq = 1000
    nonbondedcutoff = 1.0
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    # Start temp, max temp, cycles, holding_steps, steps at each temp
    anneal_parameters = [300, 700, 5, 10, 500000]
    minimized_only = None # This will change to True/False and will determine how periodic box vectors are set
    acrylate_check_freq = 100
    restrain_heavys = False
    
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
                            'Annealing_NVT': {'Time taken': None, 'Simulation time': None, 'Start temp': None, 'Target temp': None, 'Cycles': None, 'Steps at plateaus': None, 'Steps at incremental temps': None, 'Timestep': None},
                            'Basic_NPT': {'Time taken': None, 'Simulation time': None, 'Temperature': None, 'Pressure': None, 'Timestep': None},
                            'Basic_NVT': {'Time taken': None, 'Simulation time': None, 'Temperature': None, 'Timestep': None},
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
        if isinstance(self, GromacsSimulation):
            return "GRO"
        else:
            return "Please specify type of simulation by calling <AmberSimulation> OR <ANISimulation> OR <GromacsSimulation> classes"
    
    def minimize_energy(self):
        """
        Function to perform energy minimization of the system using Langevin dynamics.
        
        USAGE:
            minimized_sim = simulation_object.minimize_energy()

        Returns:
            minimized_simulation_object
        """
        min_start_time = time.time()
        
        integrator = LangevinIntegrator(self.min_temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)

        try:
            platform = Platform.getPlatformByName("CUDA")
            platform.setPropertyDefaultValue("Precision", "mixed")  # or 'single', 'double'
            print("Using CUDA platform for GPU execution.")
        except Exception:
            print("CUDA not available — falling back to CPU.")
            platform = Platform.getPlatformByName("CPU")
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator, platform)
            simulation.context.setPositions(self.amb_coordinates.positions)
            
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*unit.nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)
            simulation.context.setPositions(self.ani_coordinates)

        if BuildSimulation.type_of_simulation(self) == "GRO":
            print("Atoms in GRO:", len(self.gro_coordinates.positions))
            print("Atoms in TOP:", self.gro_topology.topology.getNumAtoms())
            system = self.gro_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            simulation = app.Simulation(self.gro_topology.topology, system, integrator, platform)
            simulation.context.setPositions(self.gro_coordinates.positions)
        
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
     
    def anneal_NVT(self, simulation, start_temp=None, max_temp=None, cycles=None, quench_rate=None, steps_per_cycle=None, filename=None, save_restart=False, restart_name=None, verbose=True):
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
        if filename is None:
            filename = "_anneal_"
        if filename is not None:
            filename = f"_{filename}_"
        if start_temp is None:
            start_temp = self.anneal_parameters[0]
        if max_temp is None:
            max_temp = self.anneal_parameters[1]
        if cycles is None:
            cycles = self.anneal_parameters[2]
        if quench_rate is None:
            quench_rate = self.anneal_parameters[3]
        if steps_per_cycle is None:
            steps_per_cycle = self.anneal_parameters[4]
         
        # Extract positional info
        state = simulation.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors() # Obtain periodc box vectors of the previous step
        
        # Set up integrator
        integrator = LangevinIntegrator(start_temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)

        try:
            platform = Platform.getPlatformByName("CUDA")
            platform.setPropertyDefaultValue("Precision", "mixed")  # or 'single', 'double'
            print("Using CUDA platform for GPU execution.")
        except Exception:
            print("CUDA not available — falling back to CPU.")
            platform = Platform.getPlatformByName("CPU")
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator, platform)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)

        if BuildSimulation.type_of_simulation(self) == "GRO":
            system = self.gro_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            simulation = app.Simulation(self.gro_topology.topology, system, integrator, platform)

        # Update the xyz of each atom
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        simulation.context.setPositions(xyz)
        
        # Total up all steps of the equilibration so the rpeorters record properly
        total_steps = steps_per_cycle*cycles
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        if self.savepdb_traj == True:
            output_pdbname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp) + ".pdb" ))
            simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
        
        # DCD trajectory
        #output_dcdname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_anneal"))
        output_dcdname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp)))
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        #output_dataname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_anneal_data"))
        output_dataname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp)))
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter)
        
        increments = int((max_temp-start_temp)/quench_rate)
        steps_per_slope = int(steps_per_cycle*0.4)
        holding_steps = int(steps_per_cycle*0.1) 
        steps_at_increment = int(steps_per_slope/increments)
        if verbose==True:
            print(f"""Annealing information:
                 - Number of heating/cooling increments: {increments}
                 - Steps per temperature in-/decrease: {steps_per_slope}
                 - Holding steps at {max_temp} k: {holding_steps}
                 - Steps at heating/cooling increment: {steps_at_increment}
                 - Total simulation time {(total_steps * self.timestep):.0f} fs
                    """)

        def cycle(start_temp, max_temp, steps_at_increment, holding_steps, increments, quench_rate):
            integrator.setTemperature(start_temp)
            simulation.step(steps_at_increment)

            for i in range(increments):
                integrator.setTemperature(start_temp + (i*quench_rate))
                simulation.step(steps_at_increment)

            integrator.setTemperature(max_temp)
            simulation.step(holding_steps)
            
            for i in range(increments):
                integrator.setTemperature(max_temp - (i*quench_rate))
                simulation.step(steps_at_increment)                

            integrator.setTemperature(start_temp)
            simulation.step(holding_steps)
                
        for i in range(cycles):
            cycle(start_temp, max_temp, steps_at_increment, holding_steps, increments, quench_rate)
            
        anneal_end_time = time.time()  
        time_taken = anneal_end_time - anneal_start_time

        self.log_info['Annealing_NVT']['Time taken'] = time_taken
        self.log_info['Annealing_NVT']['Simulation time'] = total_steps*self.timestep
        self.log_info['Annealing_NVT']['Start temp'] = start_temp
        self.log_info['Annealing_NVT']['Target temp'] = max_temp
        self.log_info['Annealing_NVT']['Cycles'] = cycles
        self.log_info['Annealing_NVT']['Steps at plateaus'] = holding_steps
        self.log_info['Annealing_NVT']['Steps at incremental temps'] = steps_at_increment
        self.log_info['Annealing_NVT']['Timestep'] = self.timestep

        # Write the final structure to pdb
        self.final_pdbname = os.path.join(self.output_dir, ("final" + filename + self.filename + ".pdb"))
        with open(self.final_pdbname, 'w') as output:
            PDBFile.writeFile(simulation.topology, state.getPositions(), output)

        if save_restart == True:
            self.save_rst(simulation, restart_name=restart_name)
            
        return(simulation, (output_dataname + ".txt"))
    
    @classmethod
    def anneal_help(cls):
        """Display help information for the anneal method."""
        print(cls.anneal.__doc__)    
    def basic_NPT(self, simulation, total_steps=None, temp=None, pressure=None, filename=None, save_restart=False, restart_name=None, verbose=True):
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
        if filename is None:
            filename = "_basic_NPT_"
        if filename is not None:
            filename = f"_{filename}_"
        if total_steps is None:
            total_steps = self.total_steps
        if temp is None:
            temp = self.temp
        if pressure is None:
            pressure = self.pressure

        if verbose==True:
            print(f"""Basic npt information
                - Total steps: {total_steps}
                - Total simulation time: {(total_steps * self.timestep):.0f} fs
                - Temperature: {temp} K
                - Pressure: {pressure} atm
                """)
            
        # Extract positional info
        state = simulation.context.getState(getPositions=True, getEnergy=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors()
        
        # Set up integrator and barostat
        barostat = MonteCarloBarostat((pressure*atmosphere), (temp*kelvin)) # Define barostat (pressure, temp)
        integrator = LangevinIntegrator(temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)

        try:
            platform = Platform.getPlatformByName("CUDA")
            platform.setPropertyDefaultValue("Precision", "mixed")  # or 'single', 'double'
            print("Using CUDA platform for GPU execution.")
        except Exception:
            print("CUDA not available — falling back to CPU.")
            platform = Platform.getPlatformByName("CPU")
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            system.addForce(barostat)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator, platform)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*unit.nanometers, constraints=app.HBonds)
            system.addForce(barostat)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)

        if BuildSimulation.type_of_simulation(self) == "GRO":
            system = self.gro_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            #system = self.restrain_heavy_atoms(system, self.gro_topology.topology, self.gro_coordinates.positions)
            #print("heavy atoms restrained")
            system.addForce(barostat)
            simulation = app.Simulation(self.gro_topology.topology, system, integrator, platform)

        # We set the box vectors with the output from the the previous simulation
        #if self.minimized_only == False:       
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        
        # Update positional info
        simulation.context.setPositions(xyz)
             
        # Set initial velocities
        simulation.context.setVelocitiesToTemperature(temp*kelvin)
        
        # Set up reporters
        if self.savepdb_traj == True:
            output_pdbname = os.path.join(self.output_dir, (self.filename +  "_" + str(pressure) + filename + str(self.timestamp) + ".pdb"))
            simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
    
        # DCD trajectory
        #output_dcdname = os.path.join(directories.systems_dir, self.filename, (self.filename +  "_" + str(pressure) + "_atm_traj.dcd"))
        output_dcdname = os.path.join(self.output_dir, (self.filename +  "_" + str(pressure) + filename + str(self.timestamp)))
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        #output_dataname = os.path.join(directories.systems_dir, self.filename, (self.filename +  "_" + str(pressure) + "_atm_data"))
        output_dataname = os.path.join(self.output_dir, (self.filename +  "_" + str(pressure) + filename + str(self.timestamp)))
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter)
        simulation.step(total_steps)       

        equili_end_time = time.time()  
        time_taken = equili_end_time - equili_start_time

        self.log_info['Basic_NPT']['Time taken'] = time_taken
        self.log_info['Basic_NPT']['Simulation time'] = total_steps * self.timestep
        self.log_info['Basic_NPT']['Temperature'] = temp
        self.log_info['Basic_NPT']['Pressure'] = pressure
        self.log_info['Basic_NPT']['Timestep'] = self.timestep

        # Write the final structure to pdb
        self.final_pdbname = os.path.join(self.output_dir, ("final_" + filename +  self.filename + "_" + str(pressure) + "_atm.pdb"))
        with open(self.final_pdbname, 'w') as output:
            PDBFile.writeFile(simulation.topology, state.getPositions(), output)

        if save_restart == True:
            self.save_rst(simulation, restart_name=restart_name)
            
        return(simulation, (output_dataname + ".txt"))
    
    @classmethod
    def basic_NPT_help(cls):
        """Display help information for the equilibrate method."""
        print(cls.basic_NPT.__doc__) 
        
    def basic_NVT(self, simulation, total_steps=None, temp=None, filename=None, save_restart=False, restart_name=None, check_acrylate=False, verbose=True):
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
        if filename is None:
            filename = "_basic_NVT_"
        if filename is not None:
            filename = f"_{filename}_"
        if total_steps is None:
            total_steps = self.total_steps
        if temp is None:
            temp = self.temp

        if verbose==True:
            print(f"""Basic npt information
                - Total steps: {total_steps}
                - Total simulation time: {(total_steps * self.timestep):.0f} fs
                - Temperature: {temp} K
                """)
            
        # Extract positional info
        state = simulation.context.getState(getPositions=True, getEnergy=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors() # Obtain periodc box vectors of the previous step
        
        # Set up integrator
        integrator = LangevinIntegrator(temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)

        try:
            platform = Platform.getPlatformByName("CUDA")
            platform.setPropertyDefaultValue("Precision", "mixed")  # or 'single', 'double'
            print("Using CUDA platform for GPU execution.")
        except Exception:
            print("CUDA not available — falling back to CPU.")
            platform = Platform.getPlatformByName("CPU")
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator, platform)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform) 

        if BuildSimulation.type_of_simulation(self) == "GRO":
            system = self.gro_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            #system = self.restrain_heavy_atoms(system, self.gro_topology.topology, self.gro_coordinates.positions)
            #print("heavy atoms restrained")
            simulation = app.Simulation(self.gro_topology.topology, system, integrator, platform)
        
        # Update positional info
        simulation.context.setPositions(xyz)
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        if self.savepdb_traj == True:
            output_pdbname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp) + ".pdb"))
            simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
        
        # DCD trajectory
        #output_dcdname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_prod_traj"))
        output_dcdname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp)))
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        #output_dataname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_prod_data"))
        output_dataname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp)))
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter) 
        if check_acrylate == False:
            simulation.step(total_steps)
        else:
            print(f"Total steps are: {total_steps}. Acrylate check steps are: {self.acrylate_check_freq}")
            number_of_checks = int(total_steps/self.acrylate_check_freq)
            print(f"Number of checks: {number_of_checks}")
            
            for i in range(number_of_checks):
                simulation.step(self.acrylate_check_freq)
                new_coords = self.save_rst(simulation, restart_name="Acrylate_check", overwrite=True)
                AcrylateReact.identify_acrylate_bonds_dist(self.topology_file, new_coords)
                # Then need to check if close enough and then update topology

        prod_end_time = time.time()
        time_taken = prod_end_time - prod_start_time
        self.log_info['Basic_NVT']['Time taken'] = time_taken
        self.log_info['Basic_NVT']['Simulation time'] = total_steps * self.timestep
        self.log_info['Basic_NVT']['Temperature'] = temp
        self.log_info['Basic_NVT']['Timestep'] = self.timestep

        # Write the final structure to pdb
        self.final_pdbname = os.path.join(self.output_dir, ("final_" + filename + self.filename + ".pdb"))
        with open(self.final_pdbname, 'w') as output:
            PDBFile.writeFile(simulation.topology, state.getPositions(), output)

        if save_restart == True:
            self.save_rst(simulation, restart_name=restart_name)
            
        return(simulation, (output_dataname + ".txt"))
        
    def thermal_ramp(self, simulation, heating=None, quench_rate=None, ensemble=None, start_temp=None, max_temp=None, total_steps=None, pressure=None, filename=None, save_restart=False, restart_name=None):
        """
        Perform a thermal ramp (heating or cooling) on the simulation system.

        This method adjusts the system's temperature incrementally, either heating or cooling it over
        a range of temperatures in defined steps. It supports both NVT and NPT ensembles, and allows
        for detailed reporting of the simulation process.

        Args:
            simulation (app.Simulation): The simulation object on which the thermal ramp is performed.
            heating (bool): If True, the system is heated; if False, the system is cooled.
            quench_rate (int): The temperature increment (in Kelvin) between successive steps of the ramp.
            ensemble (str): The ensemble used for the ramp, either "NVT" (constant volume) or "NPT" (constant pressure).
            start_temp (float, optional): The starting temperature for the ramp. Defaults to the first value in `self.anneal_parameters`.
            max_temp (float, optional): The maximum temperature for the ramp. Defaults to the second value in `self.anneal_parameters`.
            total_steps (int, optional): The total number of steps for the ramp. Defaults to `self.total_steps`.
            pressure (float, optional): The pressure to use in NPT ensemble. Defaults to `self.pressure`.

        Returns:
            tuple: The updated simulation object and the path to the output data file.

        Raises:
            ValueError: If an invalid ensemble is specified.

        Notes:
            - This method uses a Langevin integrator with the specified starting temperature.
            - The simulation system can be either AMBER (AMB) or ANI type.
            - Various simulation data and trajectories are saved, including DCD and PDB files.
            - The function updates internal logs with details of the thermal ramp.

        Example:
            simulation = app.Simulation(topology, system, integrator)
            ramp = self.thermal_ramp(simulation, heating=True, quench_rate=10, ensemble="NVT", start_temp=300, max_temp=500)
        """

        if ensemble == "NVT" or ensemble == "NPT":
            pass
        else:
            print("Please specify 'NVT' or 'NPT' ensemble for the thermal ramp")
            return()
        if heating == None:
            print("Please specifiy 'True' for heating or 'False' for cooling")
            return()
        if quench_rate == None:
            print("Please specify a quench rate as an integer")
            return()
        thermal_ramp_start_time = time.time()
        if filename is None:
            filename = "_thermal_ramp_"
        if filename is not None:
            filename = f"_{filename}_"
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
        #print("Periodic box vectors are: ", vx, vy, vz)
        
        # Set up integrator
        integrator = LangevinIntegrator(start_temp*kelvin, self.friction_coeff/picoseconds, self.timestep*femtoseconds)

        try:
            platform = Platform.getPlatformByName("CUDA")
            platform.setPropertyDefaultValue("Precision", "mixed")  # or 'single', 'double'
            print("Using CUDA platform for GPU execution.")
        except Exception:
            print("CUDA not available — falling back to CPU.")
            platform = Platform.getPlatformByName("CPU")
        
        if BuildSimulation.type_of_simulation(self) == "AMB":
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            if ensemble == "NPT":
                barostat = MonteCarloBarostat((pressure*atmosphere), ((start_temp if heating==True else max_temp)*kelvin))
                system.addForce(barostat)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator, platform)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform)

        if BuildSimulation.type_of_simulation(self) == "GRO":
            system = self.gro_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            if ensemble == "NPT":
                barostat = MonteCarloBarostat((pressure*atmosphere), ((start_temp if heating==True else max_temp)*kelvin))
                system.addForce(barostat)
            simulation = app.Simulation(self.gro_topology.topology, system, integrator, platform)

        # Update the xyz of each atom
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        simulation.context.setPositions(xyz)

        method = "heat" if heating==True else "cool"
        output_filename = self.filename + filename + method + str(self.timestamp)
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        if self.savepdb_traj == True:
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

        if heating:
            if max_temp <= start_temp:
                raise ValueError(f"Heating selected, but max_temp ({max_temp}) <= start_temp ({start_temp})")
            incremental_temps = np.arange(start_temp, max_temp + quench_rate, abs(quench_rate)).tolist()
        else:
            if start_temp <= max_temp:
                raise ValueError(f"Cooling selected, but start_temp ({start_temp}) <= max_temp ({max_temp})")
            incremental_temps = np.arange(start_temp, max_temp - abs(quench_rate), -abs(quench_rate)).tolist()

        # Sanity check
        if len(incremental_temps) == 0:
            raise ValueError(
                f"No temperature increments generated! "
                f"Check start_temp={start_temp}, max_temp={max_temp}, quench_rate={quench_rate}"
            )

        steps_at_increment = int(total_steps / len(incremental_temps))

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

        # Write the final structure to pdb
        self.final_pdbname = os.path.join(self.output_dir, ("final_" + output_filename + ".pdb"))
        with open(self.final_pdbname, 'w') as output:
            PDBFile.writeFile(simulation.topology, state.getPositions(), output)

        if save_restart == True:
            self.save_rst(simulation, restart_name=restart_name)
            
        return(simulation, (output_dataname + ".txt"))
        
    def strain(self, simulation, total_steps=None, temp=None, filename=None, save_restart=False, restart_name=None):
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
        if filename is None:
            filename = "_strain_"
        if filename is not None:
            filename = f"_{filename}_"
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
            system = self.amb_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            barostat = MonteCarloBarostat((self.pressure*atmosphere), (self.temp*kelvin))
            system.addForce(barostat)
            simulation = app.Simulation(self.amb_topology.topology, system, integrator)
        
        if BuildSimulation.type_of_simulation(self) == "ANI":
            system = self.potential.createSystem(self.ani_topology, nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=app.HBonds)
            platform = PLatform.getPlatformByName('CUDA')
            simulation = app.Simulation(self.ani_topology, system, integrator, platform) 

        if BuildSimulation.type_of_simulation(self) == "GRO":
            system = self.gro_topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=self.nonbondedcutoff*nanometers, constraints=HBonds)
            simulation = app.Simulation(self.gro_topology.topology, system, integrator)
        
        # Update positional info
        simulation.context.setPositions(xyz)
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        if self.savepdb_traj == True:
            output_pdbname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp) + ".pdb"))
            simulation.reporters.append(app.PDBReporter(output_pdbname, self.reporter_freq))
        
        # DCD trajectory
        #output_dcdname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_prod_traj"))
        output_dcdname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp)))
        dcdWriter = DcdWriter(output_dcdname, self.reporter_freq)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        #output_dataname = os.path.join(directories.systems_dir, self.filename, (self.filename + "_prod_data"))
        output_dataname = os.path.join(self.output_dir, (self.filename + filename + str(self.timestamp)))
        dataWriter = DataWriter(output_dataname, self.reporter_freq, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter) 
        
        total_strain = 0.1
        strain_steps = 100
        strain_increment = total_strain/strain_steps

        for step in range(strain_steps):
            new_vx = vx * (1 + strain_increment * (step+1))
            simulation.context.setPeriodicBoxVectors(new_vx, vy, vz)
            simulation.step(1000)

        prod_end_time = time.time()
        time_taken = prod_end_time - prod_start_time
        #self.log_info['Production']['Time taken'] = time_taken
        #self.log_info['Production']['Simulation time'] = total_steps * self.timestep
        #self.log_info['Production']['Temperature'] = temp
        #self.log_info['Production']['Timestep'] = self.timestep

        if save_restart == True:
            self.save_rst(simulation, restart_name=restart_name)
            
        return(simulation, (output_dataname + ".txt"))
        
    def save_rst(self, simulation, restart_name=None, overwrite=False):
        if restart_name == None:
            restart_name = f'final_{self.filename}'
        else:
            restart_name = f'{self.filename}_{restart_name.strip()}'
        import parmed as pmd
        from simtk import unit
        # Load original AMBER top and coord into parmed
        parm = pmd.load_file(self.topology_file, self.coordinates_file)
        # Final state
        final_state = simulation.context.getState(getPositions=True, getVelocities=True)

        parm.coordinates = final_state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
        parm.velocities = final_state.getVelocities(asNumpy=True).value_in_unit(unit.angstroms/unit.picoseconds)
    
        # save rst7 to a new folder
        restart_folder = os.path.join(self.manager.systems_dir, restart_name)
        if not os.path.exists(restart_folder):
            os.makedirs(restart_folder)
        rst_filename = os.path.join(restart_folder, f"{restart_name}.rst7")

        if overwrite == True:
            parm.save(rst_filename, format='rst7', overwrite=overwrite)
        else:
            parm.save(rst_filename, format='rst7')

        # Finally copy the topology
        new_top_name = os.path.join(restart_folder, f"{restart_name}.prmtop")
        shutil.copy(self.topology_file, new_top_name)
        return(rst_filename)
    
    @classmethod
    def production_run_help(cls):
        """Display help information for the production run method."""
        print(cls.production_run.__doc__)

    def restrain_heavy_atoms(self, system, topology, positions):
        force = openmm.CustomExternalForce("1000*(x-x0)^2 + 1000*(y-y0)^2 + 1000*(z-z0)^2")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")

        for atom in topology.atoms():
            if atom.element.symbol != "H":  # restrain heavy atoms
                pos = positions[atom.index]
                force.addParticle(atom.index, [pos.x, pos.y, pos.z])

        system.addForce(force)
        return system
        
        
    def __repr__(self):
        """
        Return a string representation of the simulation parameters for debugging and logging purposes.

        This method prints and returns a formatted string that includes key simulation parameters:
        pressure, temperature, timestep, friction coefficient, total steps, and reporter frequency.

        Returns:
            str: A formatted string representing the simulation parameters.
        """
        print("Simulation parameters: ('{}', '{}', '{}, {}, {}, {}')".format(self.pressure, self.temp, self.timestep, self.friction_coeff, self.total_steps, self.reporter_freq))
        return("Simulation parameters given in the following format: ('{}', '{}', '{}, {}, {}, {}')".format("pressure", "temperature", "timestep", "friction coefficient", "total steps", "reporter freqeuncy"))
    
    def __str__(self):
        """
        Return a user-friendly string representation of the simulation object.

        This method provides a concise description of the simulation, including its associated filename.

        Returns:
            str: A string describing the simulation object.
        """
        return 'Simulation object of - {}'.format( self.filename)
    
    @classmethod
    def display_start_time(cls):
        """
        Display the simulation start time.

        This class method prints the timestamp indicating when the simulation was initiated.

        Args:
            cls (type): The class on which the method is called.

        Returns:
            None
        """
        print("Simulation initiated at: ", cls.timestamp)

    @classmethod
    def savepdb_trajectories(cls, boolean):
        """
        Set whether to save simulation trajectories in PDB format.

        This method allows the user to toggle the saving of simulation trajectories
        in both PDB and DCD formats. If enabled, trajectories are saved in both formats;
        otherwise, only DCD format is used.

        Args:
            boolean (bool): A flag indicating whether to save PDB trajectories.
                        - `True` to save both PDB and DCD formats.
                        - `False` to save only DCD format.

        Raises:
            ValueError: If the provided argument is not a boolean.

        Returns:
            None

        Example:
            # Enable PDB trajectory saving
            Simulation.savepdb_trajectories(True)

            # Disable PDB trajectory saving
            Simulation.savepdb_trajectories(False)
        """
        if not isinstance(boolean, bool):
            print("Please pass True or False to this function.")
            return
        cls.savepdb_traj = boolean
        if boolean:
            print("Simulation will save trajectories in both .pdb and .dcd format.")
        else:
            print("Simulation will save trajectories in .dcd format only.")
    
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
        """
        Set the non-bonded cutoff distance for the simulation.

        This method allows the user to update the non-bonded cutoff distance, which
        determines how far the simulation considers interactions between non-bonded
        particles.

        Args:
            nonbondedcutoff (float): The new cutoff distance in nanometers.

        Raises:
            ValueError: If the provided argument is not a float.

        Returns:
            None

        Example:
            Simulation.set_nonbondedcutoff(5.0)
        """
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
            format_str = "Expected format: [start_temp, max_temp, cycles, quench_rate, steps_per_cycle]"
            raise ValueError(f"Invalid parameters provided. {format_str}")
        else: 
            # new_anneal_parameters = [start_temp, max_temp, cycles, holding_steps]
            cls.anneal_parameters = new_anneal_parameters
            print("Anneal parameters set.")
            print("Starting temperature is: ", str(new_anneal_parameters[0]))
            print("Target temperature is: ", str(new_anneal_parameters[1]))
            print("Number of annealing cycles is: ", str(new_anneal_parameters[2]))
            print("The quench rate is: ", str(new_anneal_parameters[3]))
            print("The number of steps per cycle is: ", str(new_anneal_parameters[4]))

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
        png_file_name = data_file.rsplit(".", 1)[0] + ".png"
        
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

class GromacsSimulation(BuildSimulation):
    """
    A class representing an AMBER molecular dynamics simulation.

    This class is designed to initialize and manage an AMBER simulation using specified 
    topology and coordinates files. It inherits from the BuildSimulation class.

    Attributes:
        manager (str): A manager object that handles directories and file paths for simulations.
        filename (str): The base name of the topology file without the extension.
        amb_coordinates (app.AmberInpcrdFile): An object representing the AMBER coordinates file.
        amb_topology (app.AmberPrmtopFile): An object representing the AMBER topology file, 
                                            which also includes periodic box vectors.

    Methods:
        __new__(cls, *args, **kwargs):
            Creates a new instance of AmberSimulation. Ensures the correct number of arguments 
            are provided before creating an instance.
        __init__(self, manager, topology_file, coordinates_file):
            Initializes an AmberSimulation object with the specified manager, topology file, 
            and coordinates file.
        __str__(self):
            Returns a user-friendly string representation of the AmberSimulation object.

    Raises:
        TypeError: If an incorrect number of arguments are provided to the __new__ method.

    Example:
        # Create a new AmberSimulation object
        sim = AmberSimulation(manager, 'path/to/topology.prmtop', 'path/to/coordinates.inpcrd')

        # Print the string representation
        print(sim)
    """
    def __new__(cls, *args, **kwargs):
        if len(args) != 3:
            raise TypeError(
                f"GromacsSimulation expected 3 arguments, but {len(args)} were given. "
                "Usage: 'sim = GromacsSimulation(manager, topology_file, coordinates_file)'"
            )
        return super().__new__(cls)

    def __init__(self, manager, topology_file, coordinates_file):
        self.manager = manager
        self.filename = os.path.basename(topology_file).split('.')[0]

        super().__init__(manager, self.filename)
        self.coordinates_file = coordinates_file
        self.topology_file = topology_file

        # Load with ParmEd
        self.gro_coordinates = GromacsGroFile(coordinates_file)
        self.gro_topology = GromacsTopFile(topology_file, periodicBoxVectors=self.gro_coordinates.getPeriodicBoxVectors())

class AmberSimulation(BuildSimulation):
    """
    A class representing an AMBER molecular dynamics simulation.

    This class is designed to initialize and manage an AMBER simulation using specified 
    topology and coordinates files. It inherits from the BuildSimulation class.

    Attributes:
        manager (str): A manager object that handles directories and file paths for simulations.
        filename (str): The base name of the topology file without the extension.
        amb_coordinates (app.AmberInpcrdFile): An object representing the AMBER coordinates file.
        amb_topology (app.AmberPrmtopFile): An object representing the AMBER topology file, 
                                            which also includes periodic box vectors.

    Methods:
        __new__(cls, *args, **kwargs):
            Creates a new instance of AmberSimulation. Ensures the correct number of arguments 
            are provided before creating an instance.
        __init__(self, manager, topology_file, coordinates_file):
            Initializes an AmberSimulation object with the specified manager, topology file, 
            and coordinates file.
        __str__(self):
            Returns a user-friendly string representation of the AmberSimulation object.

    Raises:
        TypeError: If an incorrect number of arguments are provided to the __new__ method.

    Example:
        # Create a new AmberSimulation object
        sim = AmberSimulation(manager, 'path/to/topology.prmtop', 'path/to/coordinates.inpcrd')

        # Print the string representation
        print(sim)
    """
    def __new__(cls, *args, **kwargs):
        # Check number of arguments provided
        # The __new__ method essentially creates a dummy instance before actually creating the real instance 
        if len(args) != 3:
            raise TypeError(
                f"AmberSimulation expected 3 arguments, but {len(args)} were given. "
                "Usage: 'sim = AmberSimulation(manager, topology_file, coordinates_file)'"
            )
        # Call __new__ to continue object creation if the argument count is correct
        return super().__new__(cls)
        
    def __init__(self, manager, topology_file, coordinates_file): 
        self.manager = manager
        self.filename = os.path.basename(topology_file).split('.')[0]
        
        # Inherit attributes defined in the parent class and pass directories to the BuildSimulation constructor
        super().__init__(manager, self.filename)
        self.coordinates_file = coordinates_file
        self.topology_file = topology_file
        self.amb_coordinates = app.AmberInpcrdFile(coordinates_file)
        self.amb_topology = app.AmberPrmtopFile(topology_file, periodicBoxVectors=self.amb_coordinates.boxVectors)

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
        """
        Set the potential function for the simulation.

        This method allows the user to change the potential function used in the simulation. 
        It can be used to swap between different potential models, such as switching to 
        an ANI (Artificial Neural Network) potential.

        Args:
            potential (object): The new potential function or object to be used in the simulation.

        Returns:
            None

        Example:
            # Set a new potential function
            AmberSimulation.set_potential(new_potential)

        Note:
            This method currently serves as a placeholder for switching to ANI 1 potential, 
            but it can be extended for other potential models as needed.
        """
        cls.potential = potential
        print("Potential set to: ", str(potential))

class AcrylateReact():
    def __init__(self):
        pass

    def identify_acrylate_bonds(topology=None, coordinates=None, step_check_freq=None):
        parm = pm.load_file(topology, coordinates)

        double_bonds = []
        sp2_types = {'c2', 'ce'}

        for bond in parm.bonds:
            a1, a2 = bond.atom1, bond.atom2
            if a1.atomic_number == 6 and a2.atomic_number == 6:
                if a1.type in sp2_types and a2.type in sp2_types:
                    double_bonds.append((a1, a2))

        print(double_bonds)

    def identify_acrylate_bonds_dist(topology=None, coordinates=None, step_check_freq=None, bond_center_cutoff=6):
        parm = pm.load_file(topology, coordinates)
        sp2_types = {'c2', 'ce'}

        double_bonds = []
        bond_centers = []

        for bond in parm.bonds:
            a1, a2 = bond.atom1, bond.atom2
            if a1.atomic_number == 6 and a2.atomic_number == 6:
                if a1.type in sp2_types and a2.type in sp2_types:
                    double_bonds.append((a1, a2))
                    # Compute center of bond
                    pos1 = parm.coordinates[a1.idx]
                    pos2 = parm.coordinates[a2.idx]
                    center = 0.5 * (pos1 + pos2)
                    bond_centers.append(center)

        print(f"Found {len(double_bonds)} C=C candidate bonds.")

        # Report pairwise bond center distances
        for (i, j) in itertools.combinations(range(len(bond_centers)), 2):
            center1 = bond_centers[i]
            center2 = bond_centers[j]
            dist = np.linalg.norm(center1 - center2)
            print(f"Center distance between bond {i} and bond {j}: {dist:.3f} Å")

            if dist < bond_center_cutoff:
                print(f"  -> Close bond centers detected (distance < {bond_center_cutoff} Å)")

    def react_acrylate_groups():
        pass
        