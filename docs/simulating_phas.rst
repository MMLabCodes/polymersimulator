Building and Simulating Monodisperse PHAs: walkthrough
======================================================

.. important::
   All notebooks and code should be run from the home directory: **~/polymersimulator**.  
   Running code from other directories may cause issues with file paths and prevent required Python modules from being loaded.

This guide will cover how to build polymers and run simulations of **polyhydroxyalkanoates (PHAs)** using pre-generated parameters.

The associated notebook can be found in the main PolymerSimulator directory and is called **simulating_monodisperse_PHAs.ipynb**.

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




2: Load Python Modules
----------------------

To begin the workflow, a few modules need to be imported first:

.. code-block:: python

   from modules.sw_directories import *
   from modules.sw_build_systems import *
   import os

**sw_directories**
   A filepath manager that can load and save different parameters for systems.

**sw_build_systems**
   A module containing classes to build polymers, create systems, and save parameters.

**os**
   Provides access to file paths and the base directory to initialize the filepath manager.


3: Initialise Manager and Builder Objects
-----------------------------------------

Now that the modules are loaded, two different classes — **PolySimManage** and **BuildAmberSystems** — are used to create the **manager** and **builder** objects:

.. code-block:: python

   manager = PolySimManage(os.getcwd())
   builder = BuildAmberSystems(manager)

As mentioned, the **manager** is a filepath manager and the **builder** contains the methods for building polymers and preparing systems for simulations.

.. note::
   In Python, these objects are called **classes**.  
   Classes can contain **functions**, which are called **methods**, that define the behaviors of the objects (i.e. what functions do they carry out, how is the package structured).

4: Build a Polymer
------------------

With the modules loaded and the **manager** and **builder** objects initialized, a polymer can be built.

PHAs listed at the beginning of this tutorial have already been parameterized using AmberTools .

.. note::
   All PHAs have been parameterized with **GAFF2** [#f1]_ and **abcg2** [#f2]_ charges.  

The parameterization process at a glance:

1. Build a trimer  
2. Parameterize the trimer  
3. Create **head**, **mainchain**, and **tail** units for the trimer  
4. Save these units in files so polymers can be built 

.. note::
   Parameterizing polymers with the polymersimulator is explained in more detail in other parts of this documentation. Please refer to it for projects where the parameters for different polymers do not already exist.

To build a polymer, two things are required:

**name of the base trimer**
   For any given polymer, this is: {prefix}_trimer ; where the prefix is the name of the polymer (i.e. 3HB)

.. code-block:: python
   Examples: "4HB_trimer", "3HB_trimer", "3HHp_trimer"
   
**The desired length of the final polymer**
   The number of monomers required in the final polymer (i.e 10)

Assign these variables in Python:

.. code-block:: python

   polymer_base_name = "3HB_trimer"
   number_of_units = 10

Pass these variables to the **gen_polymer_pdb_and_params** method of the builder object and assign the output to a variable called **polymer**:

.. code-block:: python

   polymer = builder.gen_polymer_pdb_and_params(
       base_name=polymer_base_name,
       number_of_units=number_of_units
   )

.. note:: 
   If you are familiar with Python, you will notice you can pass the **polymer_base_name** and **number_of_units** directly to the function. They are defined separately here for clarity.

5: Outputs from Building a Polymer
----------------------------------

The variable **polymer** contains three outputs:

1. PDB filepath of the built polymer  
2. Amber topology filepath of the built polymer  
3. Amber coordinate filepath of the built polymer

These outputs can be examined as follows and will provide filepaths to the generated files:

.. code-block:: python

   print(f"""
      Polymer built using units parameterized for: {polymer_base_name}

      The PDB file can be found at: {polymer[0]}
      The Amber topology file can be found at: {polymer[1]}
      The Amber coordinate file can be found at: {polymer[2]}""")

For the example of 3HB, the final polymer can be visualized in VMD from the PDB file and should look similar to this:

.. image:: images/3HB_decamer.PNG

.. note::
   These new files for the contstructed polymer will be in their own folder:  
   **~polymersimulator/pdb_files/systems/3HB_10_polymer**  

   The files for building the polymer units can be found at: **~polymersimulator/pdb_files/molecules/3HB_trimer**

5.1: Polymer Naming Conventions
-------------------------------

All PHAs are parameterized from trimers (e.g., `"3HB_trimer"`, `"4HB_trimer"`, etc.).  

Built polymers follow the naming pattern:

.. code-block:: none

   {prefix}_{number_of_units}_polymer

For example, using a 3HB trimer with 10 units, the polymer name will be:

.. code-block:: none

   3HB_10_polymer

5.2: Loading Polymer Files
--------------------------

While not critical to this guide, it is useful to understand how to load individual polymer files.  

Continuing with **"3HB_10_polymer"**, the  pdb files and amber topology/coordinate files can be loaded using the **manager** object:

.. code-block:: python

   polymer_name = "3HB_10_polymer"
   pdb = manager.load_pdb_filepath(polymer_name)
   amb_top, amb_coord = manager.load_amber_filepaths(polymer_name)

These variables can then be examined similarly to before:

.. code-block:: python

   print(f"""
      Polymer built using units parameterized for: {polymer_name}

      The PDB file can be found at: {pdb}
      The Amber topology file can be found at: {amb_top}
      The Amber coordinate file can be found at: {amb_coord}""")

This method works for all polymers, molecules, and systems. Only the name needs to be passed to the appropriate method in the filepath manager.

6: Building Amorphous Systems of Polymers
-----------------------------------------

The next step is to build an amorphous system of polymers with **Polyply** [#f3]_.  
There is one issue with the current setup: the polymers were parameterized and built with **AmberTools** [#f4]_, but **Polyply** was developed to be used with **GROMACS** [#f5]_. This means the current topologies are in the wrong format and a conversion to GROMACS file formats is required.

6.1: Converting Amber Topologies to GROMACS
-------------------------------------------

AmberTools has a module called **acpype** [#f6]_ which can convert topologies from Amber → GROMACS format. 

.. note::
   A more detailed explanation of Amber → GROMACS conversion will be added to the in-depth documentation.  
   A function has been implemented in PolymerSimulator for running this conversion esaily, which is what is demonstrated in this quickstart guide.

The function only requires inputs that have already been defined:

- Polymer name  
- Polymer topology  
- Polymer coordinates

This conversion is carried out with:

.. code-block:: python

   builder.run_acypype(name=polymer_name, top=amb_top, coord=amb_coord)

6.2: Building a System with Polyply
-----------------------------------

.. note::
   A more detailed explanation of this function will be added to the documentation.  
   For this quickstart guide, only the usage of the function is demonstrated.

Once the polymer has been converted to GROMACS format, multiple instances of this polymer can be packed using Polyply.  
A function called **run_polyply** within the builder object performs this task.  

The arguments required are a list of polymer names and a corresponding list of the number of each polymer.  

For example, to pack a system of **25 3HB_10_polymers**, use:

.. code-block:: python

   polymer_names = ["3HB_10_polymer"]
   number_of_polymers = [25]

These can be passed to **builder.run_polyply**:

.. code-block:: python

   system_name, gro_top, gro_coord, gro_itp = builder.run_polyply(
       polymer_names=polymer_names,
       num_poly=number_of_polymers
   )

The final system, visualised in vmd, will look similar to this (each colour is corresponds to either a head, mainchain or tail unit):

.. image:: images/3HB_10_poylymer_10_amorph.PNG

There are some noticeable *floating* atoms and bonds, this is nothing to worry about and these are atoms and bonds that lie accross the periodic boudnary conditions.


6.3 Issues with polyply starting systems
----------------------------------------

When running a simulation with a system generated with polyply, a common error is encountered:

.. code-block:: python

   OpenMMException: Particle coordinate is NaN.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan

This error typically occurs because atoms are too close to eachother and create massive repulsive forces (leading to a near infinite term for these forces in the Lennard-Jones potential). This is an artifact from the polyply packing process. The polymers are packed coarsely as minimized representations before being expanded during energy minimization. 

.. list-table::
   :widths: 50 50
   :header-rows: 0

   * - .. image:: images/polyply_out.PNG
          :width: 90%
     - .. image:: images/polyply_em.PNG
          :width: 90%

In the above images the **left** shows the polyply output and the condensed polymers. On the **right**, after energy minimization, this is what the final polymer system looks like. Due to this minimized representation --> packing --> expansion approach, there can be some unwanted steric clashes within the system leading to a system that produces a NaN error. This effect is non-existent at very low densities but quickly becomes an issue when trying to pack high-density systems. With this in mind, a default target of a system with 0.75 g/ml density is given as the desired structure.

To avoid this, a series of extremely short simulations can be carried wtih the **builder.find_polyply_starting_structrue** method. The outputs are the same as **builder.run_polyply** to avoid any confusion - except the generated files have been succesfully used to run an openmm simulation.

.. code-block:: python

   find_polyply_starting_structure(polymer_names=polymer_names, 
      num_poly=number_of_polymers, 
      max_attempts=100)

.. note::
   Test this in your own pc with a very small system. Issues should only be encoutnerred with extremely large systems but it is always worth it to be sure.

7. Running simulations
----------------------
Now a starting structure has been found, simulations can now be ran.

Import simulation module and intialise simulation
-------------------------------------------------

The first step is to import the **sw_openmm** module. This contains all of the methods to run simulations, set parameters and visualise output data.

.. code-block:: python

   from modules.sw_openmm import *

Then the toology and coordinate files can be loaded using the name of the system and the manager object.

.. code-block:: python

   system_name = "3HB_10_polymer_25_amorph"
   gro_top, gro_coord = manager.load_gromacs_filepaths(system_name)

This returns the two files required to intialise the openmm simulation - **gro_top** (topology) and **gro_coord** (coordinates).

.. code-block:: python

   sim = GromacsSimulation(manager, gro_top, gro_coord)

This returns a simulation object constructed with the given files and it is ready to do some MD!

What is the simulation object?
------------------------------

Now the simulation has been defined, a lot of information is now contained in a single python object. A list of some key attributes is given below and what they contain is given below.

.. list-table:: Simulation Attributes
   :header-rows: 1
   :widths: 30 70

   * - **Attribute**
     - **Description**
   * - ``sim.filename``
     - The name of the system.
   * - ``sim.timestep``
     - The timestep of the system (fs).
   * - ``sim.temp``
     - The temperature (K).
   * - ``sim.pressure``
     - The pressure (atm).
   * - ``sim.total_steps``
     - The total number of steps in the simulation.
   * - ``sim.nonbondedcutoff``
     - The nonbonded interaction cutoff (Å).
   * - ``sim.topology_file``
     - Path to the topology file.
   * - ``sim.coordinates_file``
     - Path to the coordinates file.

In the case of parameters like **pressure**, **temperature**, **timestep**, etc.. these are default parameters and can be changed with a single line in between simulation steps.

How to change parameters?
-------------------------
.. note::
   Many of the parameters can be passed directly to function, overiding the default parameters, but the default parameters can also be changed.

.. list-table:: Functions for Changing Simulation Parameters
   :header-rows: 1
   :widths: 35 45 20

   * - **Function**
     - **Effect**
     - **Input Format**
   * - ``sim.set_total_steps(total_steps)``
     - Updates the total number of simulation steps.
     - ``int``
   * - ``sim.set_temperature(temperature)``
     - Sets the simulation temperature.
     - ``float`` (Kelvin)
   * - ``sim.set_pressure(pressure)``
     - Sets the simulation pressure.
     - ``float`` (atm)
   * - ``sim.set_timestep(timestep)``
     - Updates the timestep of the simulation.
     - ``float`` (fs)

Where are the simulation files saved?
-------------------------------------

The last thing (*i promise!*) to be aware of before running simulations is where the simulation outputs saved?

The manager object handles the creation of a simulation directory and saving the outputs to this directory and the constructed filepaths to the output directory will follow this general form.

```python

 sim_dir = f"~polymersimulator/pdb_files/systems/{system_name}/{date_and_timestamp}"   

An an example for the simulation that is being carried in this guide will look like this:

```python

   sim_dir = "~polymersimulator/pdb_files/systems/3HB_10_polymer_25_amorph/2025-01-01_0000"

.. note::
   Each simulation is given a unique timestamp so that multiple instances of the same simulation can be ran without the files ever overwriting eachother.

The path for the simulation output directory can also be printed out with:

```python

   print(sim.output_dir)

Minimizing the energy of the system
-----------------------------------

Now the simualtions is intialized with (and assigned to the variable: **sim**), different methods can be applied to the simulations. The first step is always an energy minimization and this is carried out as so:

.. code-block:: python

   min_sim = sim.minimize_energy()

An output is assigned to this method so it can be passed to the next stage of the simulation.

.. note::
   The assignment of these output variables is critical as it allows the openmm simulation methodology to be modular as any output variable from any stage can be passed to any stage.

Next steps in a simulation
--------------------------

There are various different methods that can be applied to the system now the energy has been minimized.

.. code-block:: python

   sim.basic_NPT()
   sim.basic_NVT()
   sim.annealing_NVT()
   sim.thermal_ramp()

These methods all require an output from either **min_sim** or one of the methods shown above and other arguments that will be shown as examples in the following few sections.

The **important** thing here, that has be reiterated a few times is that these methods are modular, however, there is one key difference to the **minimize_energy** method. That method only assigned a single variable to the outputs, but for every other method, two different variables should be assigned:

- Simulation object that can be passed to another step
- Path to the data file, for quick visualisation of the results

.. note::
   The order of methods that are chosen are entirely down to the user but must make sense in the context of the project. For this example, the methods should be ordered as so, 0: Energy min., 1: NPT density equilibration, 2: NVT annealing cycle, 3: NPT heating production run.

Order of simulation stages in this example
------------------------------------------

For the system in this guide **3HB_10_polymer_25_amorph** the workflow is as follows:

- Energy minimization (already carried out at this point)
- Short NPT density equilibration (ensuring the system has reached the correct density)
- Singular NVT annealing cycle (ensure any conformationl bias inferred by packing is removed from the initial system)
- Thermal ramping NPT production run (this is the production run and intends to find thermodynamic properties of the system)

Density equilibration
---------------------

This step aims to allow the system to reach the correct density (remember systems were packed to a density of 0.75 g/ml). Before running this stage, some parameters need to be set:

.. code-block:: python

   sim.set_total_steps(10000)
   sim.set_temperature(300)
   sim.set_timestep(2.0)
   sim.set_pressure(1.0)

Setting the parameters before each step is useful to ensure simulations are running as intended.

.. note::
   Note: The total number of steps for any stage in this notebook will ideally be more and equate to a much longer simulation time than show in this notebook, however, these are just examples and longer simulations should be ran in hpc.

The variable **min_sim** can now be passed to the **basic_NPT** method and is the only argument required (all parameters have already been set).

.. code-block:: python

   npt_sim, npt_sim_data = sim.basic_NPT(min_sim)

This step assigns two variables with **npt_sim** being the important one that can be passed to the next stage in the simulation. **npt_sim_data** is the path to the data file and this filepath can be retrieved by printing the variable.

.. important::
   Whilst steps have been taken to avoid "Particle coordinate is NaN" errors, there is a possibility is persists after the openmm energy minimization. If this error occurs, rerun the simulation from the beggining. This error will be something that will be something that is avoided in a future iteration by using python exceptions.

.. code-block:: python

   print(npt_sim_data)

There is an inbuilt funciton in the simulation methodology that can produce some quick graphs from this data file:

.. code-block:: python

   sim.graph_state_data(npt_sim_data)

.. note::
   These graphs are not super important for final analysis as there are many other things that will be of interest to analyse. However, they are an easy way to ensure the simulation step was working as intended by making sure things like; temperature and density evolved throughout the simulation as expected.

Annealing
---------

The final stage before running the production run is to anneal the simulation and remove any conformational bias imposed by the initial conformation of the system.

Setting parameters for annealing follows a slightly different approach to that of the previous simulation step.

.. code-block:: python

   sim.set_anneal_parameters([start_temp, target_temp, cycles, quench_rate, total_steps])

Start temp: the temperature the annealing will start at
Target temp: the temperature the annealing will reach
Cycles: the number of annealing cycles
Quench rate: how quick the temperature will in-/de-crease
Total steps: total steps for the annealing process

For this simulation the annealing parameters that will be assigned are:

Start temp: 300
Target temp: 600
Cycles: 1
Quench rate: 10
Total steps: 100000

.. code-block:: python

   sim.set_anneal_parameters([300, 600, 1, 10, 10000])

Running the simulation is then the same as in the previous **basic_NPT** simulation, where two varaibles are assigned.

.. code-block:: python

   annealed_sim, annealed_sim_data = sim.anneal_NVT(npt_sim)

And the data can be show in a similar way.

.. code-block:: python

   sim.graph_state_data(annealed_sim_data)

Thermal ramping production
--------------------------

Now the system is at the correct density and has been annealed, the production run is finally able to be ran! (YAYYY :P)

The idea here is to heat the system from 300 K to past its experimental glass transition temperature and be able to calculate this (and other thermodynamic properties). This means that a target temperature of 600 K + is ideal.

.. note::
   These graphs are not super important for final analysis as there are many other things that will be of interest to analyse. However, they are an easy way to ensure the simulation step was working as intended by making sure things like; temperature and density evolved throughout the simulation as expected.

FROM HERE

References
----------

.. [#f1] https://doi.org/10.1021/acs.jctc.5c00038
.. [#f2] https://doi.org/10.1021/acs.jctc.8b01039
.. [#f3] https://doi.org/10.1038/s41467-021-27627-4
.. [#f4] https://doi.org/10.1021/acs.jcim.3c01153
.. [#f5] https://doi.org/10.1016/j.softx.2015.06.001
.. [#f6] https://doi.org/10.1186/1756-0500-5-367
