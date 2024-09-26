# The descriptions here are also in the README.md of the main directory "polymer_simulator"

They are additionally included here for easy finding.

## md_full_sim_v2.py ##

This script carries out a specific type of simulation aimed at yielding results from which the glass transition temperature can be found. <br>

The 'recipe' of this simulation is as follows: <br>
    1: Energy minimization <br>
    2: Ramping the temperature of system from 0 - 700 K in 50 K increments - 2.0 fs timestep (NVT) <br>
    3: Equilibration of the system at 700 K (NPT) 1.0 fs timestep <br>
    4: Cooling of the system from 700 - 300 K (NPT) 1.0 fs timestep<br> 
    
It is possible that all of the timesteps in the simulation can be 2.0 fs when simulations of long polymers are considered (there should be no explosions!). <br>
However, if a solvent is in use,  the timestep is recommended to be kept at 1.0 fs for the majority of the system and it is at the users digression to <br>
what timestep is used in the second step (maybe you want to use a large timestep up until 300 K and then move to a smaller timestep.
