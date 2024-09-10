# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 10:07:52 2024

@author: danie
"""
'''
This script equilibrates the system at 1 atm, 300k for 1 ns (npt ensemble)

This script runs a production run for 2ns at 300k for 2 ns (nvt ensemble)
'''


## Importing python modules ##
from modules.sw_directories import *
from modules.sw_openmm import *
import os as os
import sys
import time

## Define the start time of the simulation ##
start_time = time.time()

## Set up manager class  ##
main_dir = os.getcwd()
manager = SnippetSimManage(main_dir)

## Interrogate the argument passed to this script ##
if len(sys.argv) < 2:
    print("Usage: python script.py <argument>")
    print("")
    print("Where <argument> is a system name for molecular dynamics")
    sys.exit(1)

## Retrieve amber topology and coordinate file ##
top, coord = manager.load_amber_filepaths(sys.argv[1])
print(f"Amber files for: {sys.argv[1]} retrieved.")
print("")

## Set up simulation class ##
simulation = AmberSimulation(manager, top, coord)
intialisation_time = time.time() - start_time
print(f"Simulation intialized in: {intialisation_time}")
print("")

## Minimize the simulation ##
print("Minizing the energy of the system...")
min_simulation = simulation.minimize_energy()
min_time = time.time() - intialisation_time
print(f"Minisation took: {min_time}")
print("")

## Set some parameters from equilibration ## 
print("Setting parameters for equilibration.")
simulation.set_total_steps(500000) # With a 2fs timestep, this will give a 1ns equilibration
simulation.set_temperature(300)
simulation.set_pressure(1)

## Equilibrate the simulation ##
print("Equilibrating the system...")
equilibrated_sim, equilibrated_sim_data = simulation.equilibrate(min_simulation)

## Save output graphs from the equilibration ##
simulation.graph_state_data(equilibrated_sim_data)
eq_time = time.time() - min_time
print(f"System equilibrated in {eq_time}")
print("")

## Set some parameters for production run ##
print("Setting parameters for production run.")
simulation.set_total_steps(1000000) # With a 2fs timestep, this will give a 2ns production run

## Run the production run ##
print("Running the production run...")
production, production_data = simulation.production_run(equilibrated_sim)

## Save the output graphs from the production run ##
simulation.graph_state_data(production_data)
prod_time = time.time() - eq_time
print(f"Production run took {prod_time}.")

simulation.write_log_csv()



