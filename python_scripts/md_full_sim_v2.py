## Importing python modules ##
from modules.sw_directories import *
from modules.sw_openmm import *
import os as os
import sys
import time
'''
This python script is for molecular dynamics.

The aim of this script is to get thermodynamic properties of polymer melts (i.e. Td, Tg, Tm)

The steps to this script are:
1 - Energy minimization
2 - Heating of the system to 700K (this needs to be above the melting point - so can be changed between polymer systems) (NVT)
3 - Equilibration of the system at 1 atm at 700K (NPT)
4 - Cooling of the system from 700K to 300K (NPT)
'''
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
print(f"Minimisation took: {min_time}")
print("")

## Now we ramp the temperature to 700K from 0K ##
print("Ramping the temperature from 0K to 700K in NVT ensemble...")
simulation.set_timestep(2.0) # We want the system to evolve quickly here
simulation.set_total_steps(2000000) # 4 ns 
heated_sim, heated_sim_data = simulation.thermal_ramp(min_sim, True, 50, "NVT", 300, 700)
simulation.graph_state_data(heated_sim_data)
print("Temperature ramped to 700k")

## Now we equilibrate the system at 700k for 2ns ##
print("Equilibrating system at 700k in NPT ensemble...")
simulation.set_timestep(1.0)
simulation.set_total_steps(2000000)
simulation.set_temperature(700)
eq_sim, eq_sim_data = simulation.equilibrate(heated_sim)
simulation.graph_state_data(eq_sim_data)
print("System equilibrated")

## Now we cool the system slowly in NVT from 700K to 300K ##
print("Cooling the system...")
simulation.set_total_steps(10000000)
cooled_sim, cooled_sim_data = simulation.thermal_ramp(eq_sim, False, 20, "NPT", 300, 700)
simulation.graph_state_data(cooled_sim_data)
print("Temperature cooled to 300K")



