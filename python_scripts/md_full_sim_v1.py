import sys
import os

modules_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "modules_with_perm"))
sys.path.append(modules_dir)

# BuildSimulation is the parent class of ANI/Amber - Simulation
from sw_openmm import *

# sw_directories is a filepath manager and will be explained in its own guide more in depth.
from sw_directories import *

# Set up filepath manager
main_dir = os.getcwd()
main_dir = main_dir[0]
directories = PolymerSimulatorDirs(main_dir)

system_name = sys.argv[1]

# Retrieve topology and coordinates
topology_file, coordinate_file = directories.retrieve_top_crds(directories, system_name)
print("Topology and coordinate files found")
print(topology_file)
print(coordinate_file)

# Initiates simulation object with given topology and coordinate files
sim = AmberSimulation(topology_file, coordinate_file)

# Print some simulation info for the output file
print("Simulation initiated")
print(str(sim))

# Prints out the parameters defined for simulations
print(repr(sim))

# Print start time
print(sim.display_start_time())

# Minimize the energy of a system
minimized_sim = sim.minimize_energy()

# Set anneal parameters
sim.set_anneal_parameters([300, 700, 3, 1000, 10])

# Run an annealing simulation and retrieve the simulation and simulation data
annealed_sim, annealed_sim_data = sim.anneal(directories, minimized_sim)

# Generate graphs for the annealing stage of the simulation
sim.graph_state_data(annealed_sim_data)

print(repr(sim))

# Set some of the simulation parameters
sim.set_temperature(300)
sim.set_total_steps(50000)
sim.set_reporter_freq(100)
sim.set_pressure(1)

# Run an equilibration simualtion and retrieve the simulation and simulation data
equilibrated_sim, equilibrated_sim_data = sim.equilibrate(directories, annealed_sim)

# Generate graphs for the equilibration stage of the simulation
sim.graph_state_data(equilibrated_sim_data)

# Run a production simualtion and retrieve the simulation and simulation data
production_sim, production_sim_data = sim.production_run(equilibrated_sim)

# Generate graphs for the production stage of the simulation
sim.graph_state_data(production_sim_data)