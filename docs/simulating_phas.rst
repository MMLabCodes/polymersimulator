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

As mentioned, the **manager** is a filepath manager and the **builder** is the workhorse for building polymers and preparing systems for simulations.

.. note::
   In Python, these objects are called **classes**.  
   Classes can contain **functions**, which are called **methods**, that define the behaviors of the objects.

4: Build a Polymer
------------------

With the modules loaded and the **manager** and **builder** objects initialized, you can now build a polymer.

PHAs listed at the beginning of this tutorial have already been parameterized using AmberTools with the GAFF2 force field and abcg2 charges.  

The parameterization process includes:

1. Build a trimer  
2. Parameterize the trimer  
3. Create **head**, **mainchain**, and **tail** units for the trimer  
4. Save these units in files so polymers can be built 

.. note::
   This parameterization process is explained in more detail at later on in the documentation at:LINK, but this guide focuses on using pre-parameterized polymers.

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
   If you are familiar with Python, you will notice you can pass the **polymer_base_name** and **number_of_units** directly to the function. They are defined separately here for clarity.

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

Step 5.2: Loading Polymer Files
-------------------------------

While not critical to this guide, it is useful to understand how to load individual polymer files.  

Continuing with **"3HB_10_polymer"**, the files can be loaded using the **manager** object:

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

Step 6: Building Amorphous Systems of Polymers
----------------------------------------------

The next step is to build an amorphous system of polymers with **Polyply**.  
There is one issue with the current setup: the polymers were parameterized and built with **AmberTools**, but **Polyply** was developed to be used with **GROMACS**. This means the current topologies are in the wrong format.

Step 6.1: Converting Amber Topologies to GROMACS
------------------------------------------------

AmberTools has a module called **acpype** which can convert topologies from Amber → GROMACS format. 

.. note::
   A more detailed explanation of Amber → GROMACS conversion will be added to the in-depth documentation.  
   A helper function has been implemented in PolymerSimulator for repeated use, which is what is demonstrated in this quickstart guide.

The function only requires inputs that have already been defined:

- Polymer name  
- Polymer topology  
- Polymer coordinates

This conversion is carried out with:

.. code-block:: python
   builder.run_acypype(name=polymer_name, top=amb_top, coord=amb_coord)

Step 6.2: Building a System with Polyply
----------------------------------------

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

These can be passed to the function as follows:

.. code-block:: python
   system_name, system_top, system_coord, system_itp = builder.run_polyply(
       polymer_names=polymer_names,
       num_poly=number_of_polymers
   )

There are additional optional arguments, but they are not covered in this quickstart guide.  
The system will be generated with a density of 0.75 g/mL by default.