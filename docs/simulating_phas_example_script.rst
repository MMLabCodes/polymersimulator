Building and Simulating Monodisperse PHAs: Example scripts
==========================================================

This notebook contains the same code detailed here: https://polymersimulator.readthedocs.io/en/latest/simulating_phas.html . If you are not familiar with this code, it is advised that you work through the linked documentation.

In this section of the documentation, the explanations have been stripped back and instead provides examples of how this code can be built up to be used as an entire python script - with the end goal of being able to reproduce similar systems with different polymers.

A series of different things are covered here:

- Script to build polymers
- Script to build multiple polymers at once
- Script to build systems containing one type of polymer
- Script to build systems containing multiple types of polymer
- Script to run simulations
- Script to carry out all of these tasks in one

If you are unsure how to launch jupyter notebook, refer here: https://polymersimulator.readthedocs.io/en/latest/installation.html#launching-jupyter-notebook

.. note::
   PHA examples provided here are intended as a starting point and parameters are provided. You can adapt the workflow to other polymer systems using the more detailed instructions later on in the documentation.

1: Taking a Look at Pre-Parameterized PHAs
------------------------------------------

Below is a list of all PHAs that are already parameterized and available upon cloning the repository.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Image
     - Name
   * - .. image:: images/3HB_trimer.PNG
          :width: 100px
          :align: center
     - **3HB** (3-hydroxybutyrate)
   * - .. image:: images/4HB_trimer.PNG
          :width: 100px
          :align: center
     - **4HB** (4-hydroxybutyrate)
   * - .. image:: images/3HV_trimer.PNG
          :width: 100px
          :align: center
     - **3HV** (3-hydroxyvalerate)
   * - .. image:: images/6HHx_trimer.PNG
          :width: 100px
          :align: center
     - **6HHx** (6-hydroxyhexanoate)
   * - .. image:: images/7HHp_trimer.PNG
          :width: 100px
          :align: center
     - **7HHp** (7-hydroxyheptanoate)
   * - .. image:: images/2HB_trimer.PNG
          :width: 100px
          :align: center
     - **2HB** (2-hydroxybutyrate)
   * - .. image:: images/2HP_trimer.PNG
          :width: 100px
          :align: center
     - **2HP** (2-hydroxypropionate)
   * - .. image:: images/3HHx_trimer.PNG
          :width: 100px
          :align: center
     - **3HHx** (3-hydroxyhexanoate)
   * - .. image:: images/3HHp_trimer.PNG
          :width: 100px
          :align: center
     - **3HHp** (3-hydroxyheptanoate)
   * - .. image:: images/3HO_trimer.PNG
          :width: 100px
          :align: center
     - **3HO** (3-hydroxyoctanoate)
   * - .. image:: images/3HD_trimer.PNG
          :width: 100px
          :align: center
     - **3HD** (3-hydroxydecanoate)
   * - .. image:: images/3HDD_trimer.PNG
          :width: 100px
          :align: center
     - **3HDD** (3-hydroxydodecanoate)
   * - .. image:: images/3H4PhB_trimer.PNG
          :width: 100px
          :align: center
     - **3H4PhB** (3-hydroxy-4-phenylbutyrate)
   * - .. image:: images/3H5PhV_trimer.PNG
          :width: 100px
          :align: center
     - **3H5PhV** (3-hydroxy-5-phenylvalerate)
   * - .. image:: images/3H6PhHx_trimer.PNG
          :width: 100px
          :align: center
     - **3H6PhHx** (3-hydroxy-6-phenylhexanoate)
   * - .. image:: images/3H7PhHp_trimer.PNG
          :width: 100px
          :align: center
     - **3H7PhHp** (3-hydroxy-7-phenylheptanoate)
   * - .. image:: images/3H8PhO_trimer.PNG
          :width: 100px
          :align: center
     - **3H8PhO** (3-hydroxy-8-phenyloctanoate)
   * - .. image:: images/3H4mMPxPB_trimer.PNG
          :width: 100px
          :align: center
     - **3H4mMPxPB** (3-hydroxy-4-methylphenoxybutyrate)
   * - .. image:: images/3H4pMPxPB_trimer.PNG
          :width: 100px
          :align: center
     - **3H4pMPxPB** (3-hydroxy-4-methoxyphenoxybutyrate)
   * - .. image:: images/3H5BzV_trimer.PNG
          :width: 100px
          :align: center
     - **3H5BzV** (3-hydroxy-5-benzylvalerate)
   * - .. image:: images/3H5PxV_trimer.PNG
          :width: 100px
          :align: center
     - **3H5PxV** (3-hydroxy-5-phenoxyvalerate)
   * - .. image:: images/3H5pFPxV_trimer.PNG
          :width: 100px
          :align: center
     - **3H5pFPxV** (3-hydroxy-5-fluorophenoxyvalerate)
   * - .. image:: images/3H5opF2PxV_trimer.PNG
          :width: 100px
          :align: center
     - **3H5opF2PxV** (3-hydroxy-5-ortho,para-difluorophenoxyvalerate)
   * - .. image:: images/3H6mMpXHx_trimer.PNG
          :width: 100px
          :align: center
     - **3H6mMpXHx** (3-hydroxy-6-methylphenoxyhexanoate)
   * - .. image:: images/3H6pMPxHx_trimer.PNG
          :width: 100px
          :align: center
     - **3H6pMPxHx** (3-hydroxy-6-methoxyphenoxyhexanoate)
   * - .. image:: images/3H7PxHp_trimer.PNG
          :width: 100px
          :align: center
     - **3H7PxHp** (3-hydroxy-7-phenoxyheptanoate)
   * - .. image:: images/3H4MeV_trimer.PNG
          :width: 100px
          :align: center
     - **3H4MeV** (3-hydroxy-4-methylvalerate)
   * - .. image:: images/3H3PhP_trimer.PNG
          :width: 100px
          :align: center
     - **3H3PhP** (3-hydroxy-3-phenylpropionate)
   * - .. image:: images/3H5pMePxV_trimer.PNG
          :width: 100px
          :align: center
     - **3H5pMePxV** (3-hydroxy-5-methylphenoxyvalerate)

2: Import python modules and initliase python objects
-----------------------------------------------------

To begin the workflow, a few modules need to be imported, and **manager** and **builder** objects initialised.

.. code-block:: python

   from modules.sw_directories import *
   from modules.sw_build_systems import *
   from modules.sw_openmm import *
   import os

   manager = PolySimManage(os.getcwd())
   builder = BuildAmberSystems(manager)

**sw_directories**
   A filepath manager that can load and save different parameters for systems.

**sw_build_systems**
   A module containing classes to build polymers, create systems, and save parameters.

**sw_openmm**
   A module containing classes and methods to run simulations in openmm.

**os**
   Provides access to file paths and the base directory to initialize the filepath manager.

3: Script to build polymers
---------------------------

PHAs listed at the beginning of this tutorial have already been parameterized using AmberTools and the prefix of any of those polymers can be passed to this script.

To build a polymer, two things are required:

**name of the base trimer**
   For any given polymer, this is: {prefix}_trimer 
   where the prefix is the name of the polymer (i.e. 3HB)

**The desired length of the final polymer**
   The number of monomers required in the final polymer (i.e 10)

.. code-block:: python

   Base trimer examples: "4HB_trimer", "3HB_trimer", "3HHp_trimer", etc...
   Polymer length examples: 10, 20, 30, etc...

The example code below will generate a 3HB decamer.

.. code-block:: python

   # Name the trimer base name and the number units the final polymer will include
   prefix = "3HB"
   trimer_base_name = f"{prefix}_trimer"
   number_of_units = 10

   # Build the polymer using the trimer_base_name and number_of_units
   polymer = builder.gen_polymer_pdb_and_params(base_name=trimer_base_name, number_of_units=number_of_units)

This example will generate 3 files related to the polymer:

- .pdb file
- .rst7 file (amber coordinates)
- .prmtop file (amber topology)

3.1: Script to build multiple polymers at once
----------------------------------------------

The methodology can easily be altered to build multiple polymers at once by implementing it into a python **for loop**.

As an example, decamers will be built for:

- 4HB
- 3HV
- 3H4MeV

.. note::
   You can alter the prefixes to any of the pre-parameterized PHAs or add even more polymers to the list.

.. code-block:: python

   # Name the trimer base name and the number units the final polymer will include
   prefixes = ["4HB", "3HV", "3H4MeV"]
   number_of_units = 10

   for prefix in prefixes:
       # Create trimer base name
       trimer_base_name = f"{prefix}_trimer"
    
       # Build the polymer using the trimer_base_name and number_of_units
       polymer = builder.gen_polymer_pdb_and_params(base_name=trimer_base_name, number_of_units=number_of_units)

This will build 3 different decamers for the given prefixes. If you are unsure of how to call the parameters and pdb files for the generated polymers, this is explained here: https://polymersimulator.readthedocs.io/en/latest/simulating_phas.html#loading-polymer-files

.. note::
   This example is limited to buuilding decamers. However, by using i notation in a for loop, a second list of differeing polymer lengths can be defined allowing for packing different lengths of the same type of polymer. An example of this will be shown in the final script of this guide.

4: Script to build systems - 1 type of polymer
----------------------------------------------

.. important::
   This step assumes the polymer has already been built.

This step pipes together 3 functions that;

- Loads amber filepaths for a polymer
- Converts these files into the GROMACS format
- Builds a system of a given number of these polymers using polyply

This example will build the same system as seen in the walkthrough - **25 3HB decamers**. Only 2 things need to be defined:

- Name of the polymer
- Amount of the polymer 

If you are unsure how the polymer is named, please refer here: https://polymersimulator.readthedocs.io/en/latest/simulating_phas.html#polymer-naming-conventions

.. code-block:: python

   # Define list containing 1 polymer and amount of that polymer
   polymer_names = ["3HB_10_polymer"]
   number_of_polymers = [25]

   # Retrieve amber files
   amb_top, amb_coord = manager.load_amber_filepaths(polymer_names[0])

   # Convert to gromacs format
   builder.run_acpype(name=polymer_names[0], top=amb_top, coord=amb_coord)

   # Build polyply system
    system_name, gro_top, gro_coord, gro_itp = builder.find_polyply_starting_structure(polymer_names=polymer_names, num_poly=number_of_polymers, dens=750,          max_attempts=100)

This will return the required filepaths to run a simulation. However, a simulation of this system was already shown in the walkthorugh notebook so lets build a more complex system - one that conntains more than 1 type of polymer.

4: Script to build systems - multiple types of polymer
------------------------------------------------------

This step is very similar to the previous one with two exceptions:

- More polymers and the amount of those polymers are defined
- A for loop is used to prepare these polymers iteratively

For this example, the **3HV**, **4HB** and **3H4MeV** decamers will packed together.

.. code-block:: python

   # Define list containing 1 polymer and amount of that polymer
   polymer_names = ["3HV_10_polymer", "4HB_10_polymer", "3H4MeV_10_polymer"]
   number_of_polymers = [10, 10, 10]

   for polymer in polymer_names:
       # Retrieve amber files
       amb_top, amb_coord = manager.load_amber_filepaths(polymer_names[i])

       # Convert to GROMACS format
       builder.run_acpype(name=polymer_names[i], top=amb_top, coord=amb_coord)

   # Build polyply system
   system_name, gro_top, gro_coord, gro_itp = builder.find_polyply_starting_structure(polymer_names=polymer_names, num_poly=number_of_polymers, dens=750,    max_attempts=100)

The files returned here are the ones that will used for the example simulation script and **system_name** will be passed to the next stage.

6: Running a simulation
-----------------------

So far, two different systems have been generated:

- A system of 25 3HB decamers
- A system of 10 4HB decamer, 10 3HV decamers and 10 3H4MeV decamers

A test script for a simulation will be shown for the second system and the files required for simulation can be loaded easily for this with the **system_name** variable defined in the previous step.

The simulation protocol is the same as shown in the walkthrough notebook:

- Short NPT density equilibration: this is the ensure the system reaches the correct density
- A singular NVT annealing cycle: Ensure any bias is removed from the initial structure
- Thermal ramping production run in NPT: This is the final run that is intended to find the the glass transition temperatur of the polymer system

A **while loop** is also implemented here to avoid any NaN errors that sometimes result after energy minimization - something that is explained in more detail in the walkthough notebook.

.. note::
   For more information on each simulation step and how the class for running simulations works, please refer here: https://polymersimulator.readthedocs.io/en/latest/simulating_phas.html#running-simulations

.. code-block:: python

   # Load gromacs topology and coordinates
   gro_top, gro_coord = manager.load_gromacs_filepaths(system_name)

   # While loop ensures no NaN errors are hit (at least at the beggining of the simulations)
   success = False
   while not success:
       try: 
           # Intialise simulation
           sim = GromacsSimulation(manager, gro_top, gro_coord)

           # Minimize the energy in the system
           min_sim = sim.minimize_energy()

           # Set total steps (2fs timestep)
           sim.set_total_steps(10000)

           # Run a simple NPT simulation
           npt_sim, npt_sim_data = sim.basic_NPT(min_sim)

           # Update success flag
           success = True
       except Exception as e:
           # Restart the initialisation step if NaN error was encountered
           print(f"""Restarting simulation, minimized eometry imposing too many forces..
    
           The error is printed below:
    
           {e}""")

   # Visaulise data from the npt sim
   sim.graph_state_data(npt_sim_data)

   # Set annealing parameters
   sim.set_anneal_parameters([300, 600, 1, 10, 10000])

   # Anneal the simulation
   annealed_sim, annealed_sim_data = sim.anneal_NVT(npt_sim)

   # Visualise data from annealing
   sim.graph_state_data(annealed_sim_data)

   # Heat the simulation
   heated_sim, heated_sim_data = sim.thermal_ramp(annealed_sim, heating=True, quench_rate=10, ensemble="NPT", start_temp=300, max_temp=600, total_steps=10000)

   # Visaulise the data from the heating stage
   sim.graph_state_data(heated_sim_data)

Hopefully this runs without errors (it should!) but this is the final stage of going from pre-parameterized units to a final full simulation - the final step is to put everthying together.

6: Final script - Pre-parameterized polymer to a simulation
-----------------------------------------------------------

Various scripts for various stages have been shown as individual entities, the final stage is to put them all together. The idea is to go from a few simple inputs to a fully fledged simulation in one smash of the enter key. 

For this, 3 things will need to be defined:

- Prefixes of the polymers
- Number of units in each polymer
- Number of polymers in the final system

.. note::
   This example utilises an i notation for loop. This allows for different lengths of the same type of polymer to packed into the same system. It is a little bit irrelavent here as each polymer packed is a decamer, but it provides an easy route to simulating polydisperse systems.

.. code-block:: python

   from modules.sw_openmm import *
   from modules.sw_directories import *
   from modules.sw_build_systems import *
   import os as os

   manager = PolySimManage(os.getcwd())
   builder = BuildAmberSystems(manager)

   # Name the trimer base name and the number units the final polymer will include
   prefixes = ["4HB", "3HV", "3H4MeV"]
   number_of_units = [10, 10, 10]
   number_of_polymers = [10, 10, 10]

   polymer_names = []

   for i in range(len(prefixes)):
       # Create trimer base name
       trimer_base_name = f"{prefixes[i]}_trimer"
    
       # Build the polymer using the trimer_base_name and number_of_units
       polymer = builder.gen_polymer_pdb_and_params(base_name=trimer_base_name, number_of_units=number_of_units[i])

       # Create polymer names list
       polymer_names.append(f"{prefixes[i]}_{number_of_units[i]}_polymer")

   for polymer in polymer_names:
       # Retrieve amber files
       amb_top, amb_coord = manager.load_amber_filepaths(polymer_names[i])

       # Convert to GROMACS format
       builder.run_acpype(name=polymer_names[i], top=amb_top, coord=amb_coord)

   # Build polyply system
   system_name, gro_top, gro_coord, gro_itp = builder.find_polyply_starting_structure(polymer_names=polymer_names, num_poly=number_of_polymers, dens=750, max_attempts=100)

   # Load gromacs topology and coordinates
   gro_top, gro_coord = manager.load_gromacs_filepaths(system_name)

   success = False
   while not success:
       try: 
           # Intialise simulation
           sim = GromacsSimulation(manager, gro_top, gro_coord)

           # Minimize the energy in the system
           min_sim = sim.minimize_energy()

           # Set total steps (2fs timestep)
           sim.set_total_steps(10000)

           # Run a simple NPT simulation
           npt_sim, npt_sim_data = sim.basic_NPT(min_sim)

           # Update success flag
           success = True
       except Exception as e:
           print(f"""Restarting simulation, minimized eometry imposing too many forces..
    
           The error is printed below:
    
           {e}""")

   # Visaulise data from the npt sim
   sim.graph_state_data(npt_sim_data)

   # Set annealing parameters
   sim.set_anneal_parameters([300, 600, 1, 10, 10000])

   # Anneal the simulation
   annealed_sim, annealed_sim_data = sim.anneal_NVT(npt_sim)

   # Visualise data from annealing
   sim.graph_state_data(annealed_sim_data)

   # Heat the simulation
   heated_sim, heated_sim_data = sim.thermal_ramp(annealed_sim, heating=True, quench_rate=10, ensemble="NPT", start_temp=300, max_temp=600, total_steps=10000)

   # Visaulise the data from the heating stage
   sim.graph_state_data(heated_sim_data)  
