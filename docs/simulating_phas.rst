Building and Simulating PHAs
============================

This guide will cover how to build polymers and run simulations of polhydroxyalkanotes from pre-generated parameters.

The associated notebook can be found in the main polymersimulator directory and is called **simulating_polyhydroxyalkanotes.ipynb**.

.. note::
   PHA examples provided here are intended as a starting point. You can adapt the workflow to work other polymer systems. Using the more detailed instructions in the documentation

Step 1: Taking a look at PHAs already parameterized
---------------------------------------------------

Below is a list of all PHAs already parameterized and available upon cloning the repository.

TABLE OF PHAS

Step 2: Load python modules
---------------------------

To begin the workflow, a few modules need to be imported first.

.. code-block:: python
   from modules.sw_directories import *
   from modules.sw_build_systems import *
   import os as os

**sw_directories**: A filepath manager that can load and save different parameters for systems.
**sw_build_systems**: A modules containing different classes to build polymers, create systems and save parameters.
**os**: Obtains filepath to base directory to initialise filepath manager.

Step 3: Initialise manager and builder objects
----------------------------------------------

Now the modules are loaded, 2 different classes are called **PolySimManage** and **BuildAmberSystems** to create the **manager** and **builder** objects.

.. code-block:: python
   manager = PolySimManage(os.getcwd())
   builder = BuildAmberSystems(manager)

As mentioned before, the **manager** is a filepath manager and the **builder** is the workhorse of building polymers and preparing systems for simulations.
   
.. note::

   In Python, these objects are called **classes**.  
   Classes can contain **functions**, which are called **methods**, that define the behaviors of the objects.

Step 4: Build a polymer
-----------------------

Now the modules are loaded, and the **manager** and **builder** objects intialised, a polymer can be built.

PHAs in the list given at the beggining of this tutorial already have been parameterized using the GAFF2 forcefield and abcg2 charges.

These PHAs are parameterized in the following way:

1: Build a trimer
2: Parameterize the trimer
3: Create a **head**, **mainchain** and **tail** unit for the trimer
4: Create files for these units so poylmers can be built (which are stored in the **<polymer>_trimer** folder).

.. note::

   This parameterization process will be explained in more detail at LINK but this guide will only touch on the process as pre-parameterized polymers are being utilised.

To build a polymer, two things are required:

1: The name of the base trimer, so that the filepath manager can load the polymer units (i.e. "3HB_trimer")
2: A desired length of the final polymer (i.e. 10)

These variables can be assigned in python as so:

.. code-block:: python

   polymer_base_name = "3HB_trimer"
   number_of_units = 10

These can be passed to the **gen_polymer_pdb_and_params** of the builder object.

.. code-block:: python

   polymer = builder.gen_polymer_pdb_and_params(base_name=polymer_base_name, number_of_units=number_of_units)

   


