Building and Simulating Monodisperse PHAs
=========================================

This guide will cover how to build polymers and run simulations of **polyhydroxyalkanoates (PHAs)** using pre-generated parameters.

The associated notebook can be found in the main PolymerSimulator directory and is called **simulating_polyhydroxyalkanotes.ipynb**.

.. note::
   PHA examples provided here are intended as a starting point. You can adapt the workflow to other polymer systems using the more detailed instructions in the documentation.

1: Taking a Look at Pre-Parameterized PHAs
------------------------------------------

Below is a list of all PHAs that are already parameterized and available upon cloning the repository.

TABLE OF PHAS

2: Load Python Modules
----------------------

To begin the workflow, a few modules need to be imported first:

.. code-block:: python
   from modules.sw_directories import *
   from modules.sw_build_systems import *
   import os

**sw_directories**: A filepath manager that can load and save different parameters for systems.  
**sw_build_systems**: A module containing classes to build polymers, create systems, and save parameters.  
**os**: Provides access to file paths and the base directory to initialize the filepath manager.

3: Initialise Manager and Builder Objects
-----------------------------------------

Now that the modules are loaded, two different classes — **PolySimManage** and **BuildAmberSystems** — are used to create the **manager** and **builder** objects:

.. code-block:: python
   manager = PolySimManage(os.getcwd())
   builder = BuildAmberSystems(manager)

As mentioned, the **manager** is a filepath manager, and the **builder** is the workhorse for building polymers and preparing systems for simulations.

.. note::
   In Python, these objects are called **classes**.  
   Classes can contain **functions**, which are called **methods**, that define the behaviors of the objects.

4: Build a Polymer
------------------

With the modules loaded and the **manager** and **builder** objects initialized, you can now build a polymer.

PHAs listed at the beginning of this tutorial have already been parameterized using AmberTools with the GAFF2 force field and abcg2 charges.  

The pre-parameterization process includes:

1. Build a trimer  
2. Parameterize the trimer  
3. Create **head**, **mainchain**, and **tail** units for the trimer  
4. Save these units in files so polymers can be built (stored in the **<polymer>_trimer** folder)

.. note::
   This parameterization process is explained in more detail at LINK, but this guide focuses on using pre-parameterized polymers.

To build a polymer, you need:

1. The name of the base trimer (so the filepath manager can load the polymer units, e.g., `"3HB_trimer"`)  
2. The desired length of the final polymer (e.g., `10` units)

Assign these variables in Python:

.. code-block:: python

   polymer_base_name = "3HB_trimer"
   number_of_units = 10

Pass these variables to the **gen_polymer_pdb_and_params** method of the builder object:

.. code-block:: python

   polymer = builder.gen_polymer_pdb_and_params(
       base_name=polymer_base_name,
       number_of_units=number_of_units
   )

.. note:: 
   If you are familiar with Python, you can pass the **polymer_base_name** and **number_of_units** directly to the function. They are defined separately here for clarity.

5: Outputs from Building a Polymer
----------------------------------

The variable **polymer** contains three outputs:

1. PDB filepath of the built polymer  
2. Amber topology filepath of the built polymer  
3. Amber coordinate filepath of the built polymer

These outputs can be examined as follows:

.. code-block:: python

   print(f"""
      Polymer built using units parameterized for: {polymer_base_name}

      The PDB file can be found at: {polymer[0]}
      The Amber topology file can be found at: {polymer[1]}
      The Amber coordinate file can be found at: {polymer[2]}""")

For the example of 3HB, the final polymer can be visualized in VMD from the PDB file and should look similar to this:

.. image:: images/3HB_decamer.PNG

.. note::
   These new files will be in their own folder:  
   **~polymersimulator/pdb_files/systems/3HB_10_polymer**  
   The files for building the polymer units can be found at: **~polymersimulator/pdb_files/molecules.3HB_trimer**

5.1: Polymer Naming Conventions
-------------------------------

All PHAs are parameterized from trimers (e.g., `"3HB_trimer"`, `"4HB_trimer"`, etc.).  

Built polymers follow this naming pattern:

    "{prefix}_{number_of_units}_polymer"

For the current example, this results in:

    "3HB_10_polymer"

5.2: Loading Polymer Files
--------------------------

While not critical to this guide, it is useful to understand how to load individual polymer files.  

Continuing with **"3HB_10_polymer"**, the files can be loaded
