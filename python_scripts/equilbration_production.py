# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:30:13 2024

@author: danie
"""
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sw_openmm import *
import sys
import os
'''
This python script runs a molecular dynamics  qeuilibration followed by a production simulation 
    with the following steps:
    - Energy minimization
    - NPT of 1,000,000 total steps
    - NVT of 2,000,000 total steps
    
The topology file and parameter file must be provided to these scripts as command line arguments
    and instructions on executing this python script can be found in the **README.md** file.
    
1. we must define a function that loads the topology and parameter files for the simulation. 
'''
def load_files(prmtop_file, inpcrd_file):
    """
    Load Amber files and return the corresponding objects.

    Parameters:
    - prmtop_file (str): Path to the Amber prmtop file.
    - inpcrd_file (str): Path to the Amber inpcrd file.

    Returns:
    - prmtop (app.AmberPrmtopFile): Amber prmtop file object.
    - inpcrd (app.AmberInpcrdFile): Amber inpcrd file object.
    - filename (str): Name of the prmtop file (without extension).

    This function attempts to load the specified Amber prmtop and inpcrd files.
    If successful, it returns the corresponding objects and the base filename (without extension).
    If unsuccessful, it prints an error message and exits the script.
    """
    try:
        prmtop = app.AmberPrmtopFile(prmtop_file)
        inpcrd = app.AmberInpcrdFile(inpcrd_file)
        filename = os.path.basename(prmtop_file).split(".")[0]
        return prmtop, inpcrd, filename
    except Exception as e:
        print(f"Error loading files: {e}")
        sys.exit(1)
'''
2. Now we can load the files and run an openmm molecular dynamics simulation
'''
if __name__ == "__main__":
    # Check if correct number of arguments provided
    if len(sys.argv) != 3:
        print("Error: Please provide exactly two file paths (prmtop and rst7).")
        sys.exit(1)
    
    # Topology and coordinate files provided as command line arguments
    prmtop_file = sys.argv[1]
    inpcrd_file = sys.argv[2]

    # Check if files exist
    if not (os.path.isfile(prmtop_file) and os.path.isfile(inpcrd_file)):
        print("Error: One or both files do not exist.")
        sys.exit(1)

    # Load files
    prmtop, inpcrd, filename = load_files(prmtop_file, inpcrd_file)
    '''
    3. Define function for an NPT simulation
    
    Note: This npt function is only used when a system is not annealed prior to equilibration.
        If you are using a system that needs to be annealead prior to simulation, another function must be used.
    '''
    
    def initial_npt_simulation(topology_file, coordinate_file, temp, pressure, total_steps):
        """
        Set up and perform an initial NPT simulation using OpenMM.

        Parameters:
        - topology_file (str): Path to the topology (prmtop) file.
        - coordinate_file (str): Path to the coordinate (inpcrd) file.
        - temp (int): Temperature in Kelvin - as an integer.
        - pressure (int): Pressure in atmospheres - as an integer.
        - total_steps (int): Total number of simulation steps.
            
        Returns:
        - simtk.openmm.app.Simulation: OpenMM Simulation object configured for NPT simulation.
                
        This function initializes an NPT simulation using the provided topology and coordinate files.
        It creates a system with a Monte Carlo barostat, sets up a Langevin integrator, and prepares
        reporters for trajectory and data output. The resulting OpenMM Simulation object is returned
        for further use.

        Note: Make sure to run the `simulation.step(total_steps)` to execute the simulation for the specified
        number of steps after calling this function.
        """
        # Get topology and coordinates from files
        filename = os.path.basename(topology_file).split('.')[0]
        prmtop = app.AmberPrmtopFile(topology_file) # Topology file
        inpcrd = app.AmberInpcrdFile(coordinate_file) # Coordinate file
    
        # Create system and integrator
        system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds) 
        barostat = MonteCarloBarostat((pressure*atmosphere), (temp*kelvin)) # Define barostat (pressure, temp)
        system.addForce(barostat) # Add the barostat to the system
        integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds, 2.0*femtoseconds)
    
        # Create simulation
        simulation = Simulation(prmtop.topology, system, integrator)
        simulation.context.setPositions(inpcrd.positions)
        
        # Set initial velocities
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        #   visualisation of the system and coloring specific residues 
        output_pdbname = filename + "_" + str(pressure) + "_atm_traj.pdb"
        simulation.reporters.append(app.PDBReporter(output_pdbname, 1000))
    
        # DCD trajectory
        output_dcdname = filename + "_" + str(pressure) + "_atm_traj"
        dcdWriter = DcdWriter(output_dcdname, 1000)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        output_dataname = filename + "_" + str(pressure) + "_atm_data"
        dataWriter = DataWriter(output_dataname, 1000, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter)

        return(simulation)
    '''
    4. Run the NPT equilibration for 1,000,000 steps 
    '''
    # Specify number of steps for equilibration
    total_steps = 1000000
    # Initialise NPT simulation
    simulation = initial_npt_simulation(prmtop_file, inpcrd_file, 300, 1, total_steps) # (topology file, coordinate file, temp, pressure, total steps - this is specifie din the line above)
    # Run simulation 
    simulation.step(total_steps)
    '''
    5. Define function for a molecular dynamics production run
    '''
    def production_run(topology_file, simulation, temp, total_steps):
        """
        Set up and perform a production run simulation using OpenMM.

        Parameters:
        - topology_file (str): Path to the topology (prmtop) file.
        - simulation (simtk.openmm.app.Simulation): OpenMM Simulation object with initial state.
        - temp (float): Temperature in Kelvin - as an integer.
        - total_steps (int): Total number of production run steps.

        Returns:
            - simtk.openmm.app.Simulation: OpenMM Simulation object configured for production run.

        This function continues a simulation from a given state and sets up additional reporters for
        trajectory and data output. The resulting OpenMM Simulation object is returned for further use.

        Note: Make sure to run the `simulation.step(total_steps)` to execute the production run for the specified
        number of steps after calling this function.
        """
        # Retrieve information from the simulation
        state = simulation.context.getState(getPositions=True, getEnergy=True) # Define state object 
        xyz = state.getPositions() # Obtain positions of the particles from previous step
        vx, vy, vz = state.getPeriodicBoxVectors() # Obtain periodc box vectors of the previous step
        
        # Get topology from file
        filename = os.path.basename(topology_file).split('.')[0]
        prmtop = app.AmberPrmtopFile(topology_file) # Topology file
        
        # Create system and integrator 
        system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds, 2.0*femtoseconds)
        
        # Create simulation
        simulation = Simulation(prmtop.topology, system, integrator)
        simulation.context.setPositions(xyz)
        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        
        # Set up reporters
        # PDB trajectory - this is slighlty redundant with the addition of the DCD trajectory, but it is still useful for 
        #   visualisation of the system and coloring specific residues 
        output_pdbname = filename + "_prod_traj.pdb"
        simulation.reporters.append(app.PDBReporter(output_pdbname, 1000))
        
        # DCD trajectory
        output_dcdname = filename + "_prod_traj"
        dcdWriter = DcdWriter(output_dcdname, 1000)
        simulation.reporters.append(dcdWriter.dcdReporter)
    
        # Datawriter - This is a more complete data writer than previously used, the file generated is a comma delimited text file
        output_dataname = filename + "_prod_data"
        dataWriter = DataWriter(output_dataname, 1000, total_steps)
        simulation.reporters.append(dataWriter.stateDataReporter)
    
        return(simulation)
    '''
    6. Run the production simulation
    '''
    # Specify the number of steps for production simulation
    total_steps = 2000000 
    # Intialise the simulation
    simulation = production_run(prmtop_file, simulation, 300, total_steps) # (topology file, coordinate file, temp, pressure, total steps - this is specifie din the line above)
    # Run the simulation
    simulation.step(total_steps)
    