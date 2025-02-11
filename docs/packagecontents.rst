Overview of modules
===================

With the amount of code in this package, it is useful to introduce the different python modules and functionalites within the package before diving in depth into how to use them with further documentation and tutorials.


*Note: All modules are prefixed by "sw_" in homage to the use of supercomputing wales in its initial creation.*

sw_directories
--------------

**sw_directories** is a filepath manager and is the workhorse of this package. This package is designed so that all of the code 
functions by being executed within the base directory.

In this vein, there are a range of functions that locate specific types of and/or specific files. The module is class based and the initial class instance
contains filepaths to varying directory locations around the package (which enables the file finding) so generated files can be found in an intuitive location.

There are other functionalities of this python module that will recieve documentation at a later date but this includes a class for:

- managing a DFT job within supercomputing wales
- organising files for generation of data-driven representative models for biocrude oils.


sw_basic_functions
------------------

**sw_basic_functions** includes basic functions that are seemingly random but have uses across all modules within this package.

This module simply imports other functions such as:

- estimating volume of a molecule (from SMILES//pdb//RDkit mol object)
- identifying structural information within a molecule

sw_build_systems
----------------

This is a class based module with a lot of content. It's main functionality is to build systems for molecular dynamics simulation.

There are two classes here; **BuildSystems** and **BuildAmberSystems** 

**BuildSystems** contains functions such as:

- SMILES to pdb file
- calling packmol
- retrieving xyz distances of a molecule
- assigning generated molecules unique residue codes and storing them in a database

**BuildAmberSystems** inherits all of the functionality of **BuildSystems** but uses AmberTools to carry out a range of other functionalities.

For example:

- Parameterize molecules
- Parameterize polymers
- Build polymers
- Build  systems of multiple polymers

sw_openmm
---------

This is a class based module and carries out molecular dynamics simulations. It is built on top of the existing openmm python module.
There are some legacy functions that will not be supported in this documentation. The main class in this module is the **BuildSimulation** class.

Whilst **BuildSimulation** is the main class, this class is a base for the ones you will need to use, these are **AmberSimulation** and **AniSimulation**.
These 2 classes inherit class methods to carry out molecular dynamic simulations from the **BuildSimulation** class.

**AmberSimulation** is a class that uses amber topology and coordinae files to carry out molecular dynamic simulations.

**ANISimulation** is a class that uses a pdb file to build and carry out molecular dynamic simulations - it is currently not supported//working.

**BuildSimulation** contains many different molecular dynamic functions such as:

- Energy minimization
- Simple NVT//NPT simulations
- Annealing in NVT
- Thermal ramping//deramping in NVT//NPT
- Easy ways to change simulation parameters
- Saving of output file

*Note: the **sw_directories** module handles the location of any output files*

sw_analysis
-----------

This is a class based module and is the most complex in its usage - it is best learnt with the different tutorials.

It analyses output trajectories from different stages of a simulation and has a variety of uses:

- Calculating free volume
- Calculating radius of gyration
- Calculating end to end distances of polymers
- Calculating persistence length of polymers
- Extracting coordinates of molecules for further study

sw_orca
-------

This is a classed based module and contains one class **orca_molecule**.

This class extracts information from a .csv containing DFT information and allows this information to be passed to other modules.



