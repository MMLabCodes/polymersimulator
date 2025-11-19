Analysing Polymer Simulations: Walkthrough
==========================================

.. important::
   All notebooks and code should be run from the home directory: **~/polymersimulator**.  
   Running code from other directories may cause issues with file paths and prevent required Python modules from being loaded.

This guide will cover how to set up and run some pre-defined analysis techniques on a polymer system after an MD simulation. It is also recommended that
you are somewhat familiar with MDAnalysis before utilising this code as it will help you understand any errors and the processes carried out in this guide.

The associated notebook can be found in the main PolymerSimulator directory and is called **Analysing_polymer_simulations_Walkthrough.ipynb**.

If you are unsure how to launch jupyter notebook, refer here: https://polymersimulator.readthedocs.io/en/latest/installation.html#launching-jupyter-notebook

.. note::
   These analysis methods are only valid for monodisperse homopolymer systems. Support for other types of systems that mix polymer types and lengths will be added in due course.

1. Setting up the analysis class
--------------------------------

To begin the workflow, a few modules need to be imported first:

.. code-block:: python

   from modules.sw_directories import *
   from modules.sw_analysis import *
   import os

**sw_directories**
   A filepath manager that can load and save different parameters for systems.

**sw_analysis**
   A module containing classes to set up analysis workflows and carry out pre-defined analysis

**os**
   Provides access to file paths and the base directory to initialize the filepath manager.

2: Initialise Manager and Polymer universe object
-------------------------------------------------

Now that the modules are loaded, two different classes — **PolySimManage** and **initialise_poly_analysis** — are used to create the **manager** and **universe** objects.

.. code-block:: python

   manager = PolySimManage(os.getcwd())

Now the manager object has been set up, it is important to note what arguments are required to set up the universe object. The universe object is a wrapper of a standard MDAnalysis universe, but it contains useful attributes specific to this analysis workflow (hence the existence of this wrapper) - the many attributes and functions will be detailed further on in this guide.

1. A finished polymer simulation
2. The name tag associated to the simulation stage you desire to analyse (see example figure below)
3. The name of your system
4. The name of your polymer in the system
5. The length of the polymers in the system
6. The type of files the simulation was launched with

As an example, files of a system called **3HB_10_polymer_25_amorph** can be found in the polymersimulator repository. This system is:

- 25 3HB decamers
- Amorphous (starting structure generated with Polyply)

Within the standard polymersimulator file structure all the files related to this system are found at **~/polymersimulator/pdb_files/systems/3HB_10_polymer_25_amorph**. Inside of this folder, various things can be found:

- Topology files
- Coordinate files
- Timestapped folders for any simulations

An example is shown below (although, bear in mind, there may be files present that you do not have due to testing and other methods not yet documented).
- Highlighted in blue are the general files for the system *(i.e. topologies and coordinates). 
- Highlighted in orange is a folder for a specific simulation and contains data files, trajectories and some output graphs. 
- Highlithed in yellow is the location of this file within the polymersimulator file structure.

.. image:: images/simulation_folder.PNG

It is worth explaining some of the files inside the actual simulation folder in a bit more detail:

.. image:: images/simulation_files.PNG

There are are different types of files:

- **.dcd**: these are trajectory files
- **.png**: these are outputted graphs from each simulatiuon stage
- **.txt**: these are the raw data files from each simulations tage
- **.csv**: this file is the log file, it is depreceated and the contents can be ignored - but if this file exists, it means the simulation finished

These files all follow a similar naming convention:

- **3HB_10_polymer_25_amorph**: the name of the system
- **"NVT_annealing"**: this is the name of the simulation stage with other examples highlighted in yellow
- **2025-11-10_....**: this is the specific timestamp assigned to this simulation - this exists to allow for replica simulations without overwriting previously existing data.

Now that long-winded explanation is out of the way, the base line of code to initialise an analysis workflow looks like this:

.. code-block:: python

   universe = intialise_poly_analysis(manager, system_name, polymer_name, poly_len, sim_stage_name, sim_type, sim_index)

Each of these arguments is explained in the table below:

.. list-table:: Arguments for ``initialise_poly_analysis``
   :header-rows: 1
   :widths: 20 80

   * - **Argument**
     - **Description**
   * - ``manager``
     - Manager object. 
        EXAMPLE: pre-intialised manager object. 
        TYPE: python object
   * - ``system_name``
     - Name of the system. 
        EXAMPLE: "3HB_10_polymer_25_amorph". 
        TYPE: string
   * - ``polymer_name``
     - Name of the polymer. 
        EXAMPLE: "3HB_10_polymer". 
        TYPE: string
   * - ``poly_len``
     - Length of the polyer. 
        EXAMPLE: 10 . 
        TYPE: integer
   * - ``sim_stage_name``
     - Name of the simulation stage to be analysed. 
        EXAMPLE: "cooling_NPT_cool". 
        TYPE: string
   * - ``sim_type``
     - Type of simulation files parsed to simulation engine (i.e. amber or gromacs files). 
        EXAMPLE: "GRO". 
        TYPE: string. 
        NOTE: types currently supported are "AMB" and "GRO"
   * - ``sim_index``
     - The simulation folder you want to access. 
        EXAMPLE: 0 . 
        TYPE: integer. 
        NOTE: this is useful when you have replica universes, you can utilise the same other arugments but pass 0,1,2,... as an argument here to analyse a different instance of the same simulation.

So for the exmaple system, the section I want to analyse is the production run which contains the "cooling_NPT_cool" tag (if you want to analyse another section of the simulation, just pass a different file tag, i.e. *"NVT_annealing"*) so the arguments are as follows:

- manger = manager
- system_name = "3HB_10_polymer_25_amorph"
- polymer_name = "3HB_10_polymer"
- poly_len = 10
- sim_stage_name = "cooling_NPT_cool"
- sim_type = "GRO"
- sim_index = 0

And the line of code looks like this:

.. code-block:: python

   universe = initialise_poly_analysis(manager=manager, system_name="3HB_10_polymer_25_amorph", polymer_name="3HB_10_polymer", poly_len=10, sim_stage_name="cooling_NPT_cool", sim_type="GRO", sim_index=0)

   





   







