# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:04:03 2024

@author: 983045
"""
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import os
'''
This python script runs a simple molecular dynamics simulations with the following steps:
    - Energy minimization
    - NVT of 10,000 total steps
    
The topology file and parameter file must be provided to these scripts as command line arguments
    and instructions on executing this python script can be found in the **README.md** file.
    
1. we must define a function that loads the topology and parameter files for the simulation. 
'''
def load_files(prmtop_file, inpcrd_file):
    ### IMPORTANT: THIS FUNCTION NEEDS TO BE EXTRACTED INTO BUILD SIMULATION CLASS ##
    
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
2. with the loaded files - run an openmm molecular dynamics simulation
'''
if __name__ == "__main__":
    # Check if correct number of arguments provided
    if len(sys.argv) != 3:
        print("Error: Please provide exactly two file paths (prmtop and rst7).")
        sys.exit(1)

    prmtop_file = sys.argv[1]
    inpcrd_file = sys.argv[2]

    # Check if files exist
    if not (os.path.isfile(prmtop_file) and os.path.isfile(inpcrd_file)):
        print("Error: One or both files do not exist.")
        sys.exit(1)

    # Load files
    prmtop, inpcrd, filename = load_files(prmtop_file, inpcrd_file)
    # Create system and integrator
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)

    # Create simulation
    simulation = app.Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

    # Minimize energy
    simulation.minimizeEnergy()

    # Set initial velocities
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

    # Run simulation
    simulation.reporters.append(app.PDBReporter((filename + "_output.pdb"), 1000))
    simulation.reporters.append(app.StateDataReporter((filename + "_data.csv"), 1000, step=True, potentialEnergy=True, temperature=True))
    simulation.step(10000)
