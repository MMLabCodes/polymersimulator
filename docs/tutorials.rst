Tutorials
=========

Running tutorials
-----------------

From the command line (or ubuntu in windows) if you launch jupyter notebook from the main directory of project, you should see a range of different jupyter tutorials.
Follow these alongside the tutorials in the documentation to get started with the python package.

sw_directories tutorial
-----------------------

The first tutorial is **Tutorial_1_filepath_manager** and will explain how files are organised and functions wrapped around them throughout the package.

sw_build_systems tutorial
-------------------------

This encompasses the second and third tutorials; **Tutorial_2_build_systems** and **Tutorial_3_build_amber_systems**. These encompass building molecules from SMILES strings,
parameterizing them and generating amber topology and coordinate files for molecular dynamics simulation. The building and parameterizaiton of polymers and
generation of amber files for systems containing these polymers is also covered.

sw_openmm tutorial
------------------

The fourth tutorial (**Tutorial_4_Running_Openmm_Simulations**) uses systems generated in tutorials 2 and 3 and runs molecular dynamics simulations of them. It 
contains a **quickstart guide** but also runs through the intricacies of the module.

sw_analysis tutorial
--------------------

Coming one day.....

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   sw_directories_tutorial
   sw_build_systems_tutorial
   sw_openmm_tutorial


