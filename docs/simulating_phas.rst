Building and Simulating PHAs
============================

This guide will cover how to build polymers and run simulations of **polyhydroxyalkanoates (PHAs)** using pre-generated parameters.

The associated notebook can be found in the main PolymerSimulator directory and is called **simulating_polyhydroxyalkanotes.ipynb**.

.. note::
   PHA examples provided here are intended as a starting point. You can adapt the workflow to other polymer systems using the more detailed instructions in the documentation.

Step 1: Taking a Look at Pre-Parameterized PHAs
-----------------------------------------------

Below is a list of all PHAs that are already parameterized and available upon cloning the repository.

TABLE OF PHAS

Step 2: Load Python Modules
---------------------------

To begin the workflow, a few modules need to be imported first:

.. code-block:: python
   from modules.sw_directories import *
   from modules.sw_build_systems import *
   import os

**sw_directories**: A filepath manager that can load and save different parameters for systems.  
**sw_build_systems**: A module containing classes to build polymers, create systems, and save parameters.  
**os**: Provides access to file paths and the base directory to initialize the filepath manager.

Step 3: Initialise Manager and Builder Objects
----------------------------------------------

Now that the modules are loaded, two different classes — **PolySimManage** and **BuildAmberSystems** — are used to create the **manager** and **builder** objects:

.. code-block:: python
   manager = PolySimManage(os.getcwd())
   builder = BuildAmberSystems(manager)

As mentioned, the **manager** is a filepath manager, and the **builder** is the workhorse for building polymers and preparing systems for simulations.

.. note::
   In Python, these objects are called **classes**.  
   Classes can contain **functions**, which are called **methods**, that define the behaviors of the objects.

Step 4: Build a Polymer
-----------------------

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
   If you are familar with python, you will notice that you can pass the **polymer_base_name** and **number_of_units** directly to the function. However, they have been defined seperately here for clarity.

Step 5: Outputs from building a polymer
---------------------------------------

When the polymer was built, the variable **polymer** was assigned the output. This output variable contains 3 things:

1: PDB filepath of built polymer
2: Amber topology filepath of built polymer
3: Amber coordinate filepath of built polymer

This variable can be interrogated as shown below.

.. code-block:: python

   print(f"""
      Polymer built using units parameterized for: {polymer_base_name}

     The pdb file can be found at: {polymer[0]}

      The amber topology file can be found at: {polymer[1]}

      The amber coordinate file can be found at: {polymer[2]}""")

For the example of 3HB, the final polymer can be visualised in VMD from the PDB file and should look similar to this.

.. image:: images/3HB_decamer.PNG

.. note::
   For clarity, these new files will be in their own folder. **~polymersimulator/pdb_files/systems/3HB_10_polymer** and the files for building the polymer units can be found at **~polymersimulator/pdb_files/molecules.3HB_trimer**

Step 5.1: Polymer naming conventions
------------------------------------

It is important to understand the polymer naming convention used when polymers are built.

As previoulsy mentioned, all the PHAs are parameterized from the trimers (i.e. "3HB_trimer", "4HB_trimer", etc...)

The built polymers will follow this naming pattern:

    "{prefix}_{number_of_units}_polymer"

And in the current example this will give

    "3HB_10_polymer"

Step 5.2: Loading polymer files
-------------------------------

This is not critical to the current guide, but it is important to understand how you can load individual files and how other processes in the code load different files for building systems and running simulations.

Following the naming convention for the polymers, and continuing with **"3HB_10_polymer"**, the files can be loaded by simply passing the name to some different functions found in the **manager** object.

.. code-block:: python

   polymer_name = "3HB_10_polymer"
   pdb = manager.load_pdb_filepath(polymer_name)
   amb_top, amb_coord = manager.load_amber_filepaths(polymer_name)

Then these variables can be interrogated in a simalar way to before:

.. code-block:: python

   print(f"""
      Polymer built using units parameterized for: {polymer_name}

      The pdb file can be found at: {pdb}

      The amber topology file can be found at: {amb_top}

      The amber coordinate file can be found at: {amb_coord}""")

This method of loading files is the same for all polymers, molecules and systems. Only the name needs to be passed to the appropriate method found in the filepath manager.

Step 6: Building amorphous systems of polymers
----------------------------------------------

The next step is to build an amorphous system of polymers with **polyply**. There is one issue with the current set up though.
The polymers were parameterized and built with **AmberTools** and **polyply** was developed to be used with **GROMACS**. This means 
the current toplogies are in the wrong format.

Step 6.1: Converting amber topologies to GROMACS
------------------------------------------------

Ambertools has a module called **acpypye** which can convert topologies from amber --> GROMACS format. 

.. note::
   A more detailed explanation of how to convert from amber --> GROMACS will be added to the in-depth documentation. However, 
a function for use within the polymersimulator has been developed to allow for repeated use of this functionality of AmberTools and will be the only thing shown in this quickstart guide.

The function only requires inputs that have already been defined, these are:

- Polymer name
- Polymer topology 
- Polymer coordinates

This conversion is carried out with the following code (all of the variables should already be defined if all of the code in this guide has been executed)

.. code-block:: python
   builder.run_acypype(name=polymer_name,top=amb_top, coord=amb_coord)

Step 6.2: Building a system with polyply
----------------------------------------

.. note::
   A more detailed explanation of how this function works will be added to the documentation. However, for this quickstart guide, only the use of function will be covered

Now the polymer has been converted to gromacs format, multiple instances of this polymer can be packed using polyply. A function has been created called **run_polyply** within the builder object to carry out this task. 

The arguments required for this function are a list of polymer names and a list of the amount of each polymer.

For example, to pack a system of **25 3HB_10_polymers** the arguments will be:

.. code-block:: python
   poylmer_names = ["3HB_10_polymer"]
   number_of_polymers = [10]

These can be passed to the function as so:

.. code-block:: python
   system_name, system_top, system_coord, system_itp = builder.run_polyply(polymer_names=polymer_names, num_poly=number_of_polymers)

There are some additional, optional, arguments that can be passed here, but those will not be covered in this tutorial. The important thing to note here is that the system will generated with a density of 0.75 g/ml.



