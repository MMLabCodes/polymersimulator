sw_build_systems tutorial
=========================

This module is the basis for almost all of the tasks carried out by other modules in this package. It relies on the **SnippetSimManage** class in the **sw_directories** module
for a lot of file finding//directing. 

The second and third tutorials (**tutorial_2_build_systems** and **tutorial_3_solvating_small_molecules**) provides a nice introduction into how to use this module and carry out of ranging of tasks related to building molecules and generating amber files for molecular dynamics.

This module consists of 2 classes; **BuildSystems** and **BuildAmberSystems**.

Initialising an instance of BuildSystems
----------------------------------------

The **BuildSystems** class gives access to a range of standard functions. Before going through some examples lets intialise an instance of this class, to do
this, an instance of the **SnippetSimManage** class will also have be initiated as explained in the **sw_directories tutorial**.

.. code-block:: python

   from modules.sw_directories import *
   from modules.sw_build_systems import *
   import os as os

   manager = SnippetSimManage(os.getcwd())
   builder = BuildSystems(manager)

Building molecules with BuildSystems
------------------------------------

The next stage is to build a molecule with the instance of **BuildSystems**. To do this a SMILES string and a name are required, for the purposes of the tutorial lets build caffeine. It SMILES string is: **CN1C=NC2=C1C(=O)N(C(=O)N2C)C**.

.. code-block:: python

   builder.SmilesToPDB("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine")

The arguments of this function are:

.. code-block:: python

   builder.SmilesToPDB(SMILES, path_to_generated_pdb_file)   

This class method will generate a .pdb file containing a molecule associated with the SMILES string saved at the specified location. 
So in the example above, a file called **caffeine.pdb** will now exist in my home directory. This function doesn't send any files to correct places
as it is intended as a test to see if your desired molecule can be built.

The alternative way to build a molecule is the one you will want to use instead as it will create a directory for your molecule and save a .pdb file there
in addition to generating a unique residue code for it.

.. code-block ::python

   pdb_file = builder.SmilsToPDB_GenResCode("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine")

For this function the only arguments the SMILES and name of a molecule:

.. code-block:: python

   pdb_file = builder.SmilesToPDB(SMILES, molecule_name)   
   print(pdb_file)

We do not need to specify a filepath for the pdb file to be saved in the correct place as when we initialised an instance of the BuildSystems class, all of the information from the filepath manager was also passed to it.

There are a couple of useful things we can now do with the manager instance:

.. code-block:: python

   pdb_file = manager.load_pdb_filepath("caffeine")
   res_code = manager.retrieve_rescode("caffeine")
   print(f"The residue code for caffeine is {res_code} and the path to the associated pdb file is {pdb_file}")

You can also open the pdb file and see that the rescode in the pdb file matches the one printed above.

.. image:: images/caffeine.PNG

The residue code information will now be stored in the residue code csv file which is located:

.. code-block:: python

   print(manager.residue_code_csv)

Other functionality of BuildSystems
-----------------------------------

The meat of this module is contained in the next class, but before moving on to that, there is some other functionality of the **BuildSystems** class.

Finding largest distances between points in an .xyz or .pdb file
----------------------------------------------------------------

The maximum distance between points in all 3 plance can be returned with the following class method.

.. code-block:: python

   x,y,z = builder.get_xyz_dists(pdb_file)
   print(x,y,z)

This function requires a pdb filepath (or xyz filepath) so you can use one the methods shown above that return the entire path to the pdb file.
This is useful when determining the sizes of organised systems contaiing many molecules and when generating periodic box vectors for a molecular dynamics simulation.

Aligning a molecule in a plane
------------------------------

A molecule can be aligned in a plane with the following class method.

.. code-block:: python

   # Can be "X"/"Y"/"Z"
   pdb_file = builder.align_molecule(pdb_file, "Z") 

This function requries a pdb filepath and a specification of what plane you want to align the molecule in. It utilises a principal axis calcualted in mdanalysis and 
rewrites the pdb file. You will want to check your pdb file visually as the results are not foolproof and the calculated principal axis may not always be as expected.
This function will return the original path of the pdb file, so you can keep working as seamlessly as possible.

Other
-----

There will be some other class methods you may see in the documentation that you can ignore for some cases. First of all, any instance of **BuildSystems**
will have a **manager** attribute, and this is everything from the **SnippetSimManage** class and this allows buildsystems to send and retireve files, but you could
technically do anything you could do with the manager instance with the builder instance as so:

.. code-block:: python

   builder.manager.some_method_or_attribute

There is also a class method that runs packmol that is called as follows:

.. code-block:: python

   builder.run_packmol(input_filepath)

However, you will need packmol configured in your pc (or hpc facility) to use this and this emulates the typical packmol command you would normally run the shell.

There are also a series of class methods related to handling residue codes and ensuring none are overwritten or duplicated, you can ignore these as these take place in the background when building molecules with unique residue codes.


Initialising an instance of BuildAmberSystems
---------------------------------------------

The second class is **BuildAmberSystems** and inherits all of the functionality from its parent **BuildSystems** class and anything you can do with **BuildSystems** is also possible with **BuildAmberSystems**.
This class is the proliferator of systems for molecular dynamics. It is recommended to follow **Tutorial_3_build_amber_systems**.
The first thing to do is initialise an instance of the class - this is the exact method as for **BuildSystems* but using **BuildAmberSystems** instead.

.. code-block:: python

   from modules.sw_directories import *
   from modules.sw_build_systems import *
   import os as os

   manager = SnippetSimManage(os.getcwd())
   builder = BuildAmberSystems(manager)

Parameterizing a molecule
-------------------------

The next stage is to parameterize a molecule for molecular dyanmics with amber. Before actually parameterizing it charges to need to be assigned to each atom
in the molecule. This file with the charges is a '.ac' file and it is essentially a '.pdb' file with an extra notation for the charges of each atom.

.. code-block:: python

   ac_file = builder.gen_ac_file("caffeine")
   print(ac_file)

The only argument for this class method is the name of the molecule being parameterized - due to the manager object (an instance of **SnippetSimManage**) being
passed to **BuildAmberSystems** it can find the directory easily. The filepath of the .ac file is also returned, but it is not required for anything other than ensuring the file exists and has actually been generated.

The next step is to parameterize caffeine.

.. code-block:: python

   builder.parameterize_mol("caffeine")

This class method doesnt return anything but there is a way to check if the molecule was parameterized.

.. code-block:: python

   param = builder.is_mol_parameterized("caffeine")
   print(param)

If it has been succesfully parameterized **True** will be printed and **False** if it has not been. This class method could be integrated into some automatic methods
and provide an instance to catch any errors in parameterization. You will only ever need to parameterize each molecule once and **warning** for larger molecules
this process will take significantly longer.

*Note: you may be a bit confused at this point as to where the polymers are...! That will be the next, but first it is critical to understand the functionality of each module and build understanding in a useful way.*

Generating amber parameter files for a single molecule
------------------------------------------------------

To run an amber simulation a topology and coordinate is required and this needs to be built using the parameter files. A class method exists
to carry out the generation of these files automatically.

.. code-block:: python

   system_name = builder.gen_amber_params_sing_mol("caffeine")
   print(system_name)

This will generate the topology and coordinate files for a single molecule of caffeine, this step may seem unnecessary as the parameters for caffine already exist
but this stage is required to generate periodic boundary conditions for the molecule. The system name will be printed and can be used to retrieve the files for molecular dyanmics simulation.
The system name in this example is **caffeine_sing_mol** and other class methods that build systems will return appropriately named systems.

.. code-block:: python

   top, coord = manager.load_amber_filepaths(system_name)
   print(f"The topology file for {system_name} is {top}")
   print(f"The coordinate file for {system_name} is {coord}")

The system will be a single caffeine molecule.

INSERT PIC

Generating amber parameter files for a single molecule solvated in water
------------------------------------------------------------------------

To generate amber parameters for a small molecule solvated in water, the procedure is very much the same, but a different class method is called.

.. code-block:: python

   system_name = builder.gen_amber_params_sing_mol_solvated("caffeine")

The amber files can be found using the system name in the same way as before. This system will be called **caffeine_wat_solv**. There is another consideration here though.
This class method will create pariodic box where the molecule exactly fits and solvate the free space. However, you may want to add a buffer of water around your molecule.
This is done by adding another argument to the class method above.

..code-block:: python

   system_name = builder.gen_amber_params_sing_mol_solvated("caffeine", 10)

This will add a buffer of water 10 angstroms around the caffeine molecule and this system will be called **caffeine_wat_solv_10** and, again, the amber files can be retrieved as so.

..code-block:: python

   top, coord = manager.load_amber_filepaths(system_name)
   print(f"The topology file for {system_name} is {top}")
   print(f"The coordinate file for {system_name} is {coord}")

INSERT PIC

   








   
   