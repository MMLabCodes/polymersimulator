Building and Simulating PHAs
============================

This guide walks you through building **Polyhydroxyalkanoates (PHAs)** and running molecular dynamics simulations with Polymer Simulator.  
It assumes you have already completed the :doc:`installation` steps and successfully tested your environment.

.. note::
   PHA examples provided here are intended as a starting point. You can adapt the workflow to other polymer systems.

Step 1: Prepare Input Files
---------------------------

Start by creating a working directory for your PHA project:

.. code-block:: bash

   mkdir pha_project
   cd pha_project

Obtain or generate the monomer structure for your chosen PHA.  
Supported formats include ``.pdb`` and ``.mol2`` files.

.. tip::
   You can use **antechamber** to preprocess small molecules and generate force field parameters.

Step 2: Parameterize the Monomer
--------------------------------

Use ``antechamber`` to generate charges and parameters:

.. code-block:: bash

   antechamber -i monomer.pdb -fi pdb -o monomer.mol2 -fo mol2 -c bcc -s 2

This produces a ``monomer.mol2`` file with assigned charges.

.. warning::
   If you encounter charge assignment errors, check that hydrogens are present in the input structure.

Step 3: Build the Polymer Chain
-------------------------------

Use the Polymer Simulator module to build a polymer chain from the monomer:

.. code-block:: python

   from polymersimulator.builder import PolymerBuilder

   # Create a polymer of 20 repeating units
   builder = PolymerBuilder("monomer.mol2", n_units=20)
   polymer = builder.build()

   polymer.save("pha_polymer.pdb")

This generates a PDB file of the polymer chain.

Step 4: Generate Topology and Coordinates
-----------------------------------------

Load the polymer into ``tleap`` to build the system:

.. code-block:: bash

   tleap -f leaprc.gaff

In the ``tleap`` shell:

.. code-block::

   source leaprc.gaff
   mol = loadpdb pha_polymer.pdb
   saveamberparm mol pha.prmtop pha.inpcrd
   quit

You now have the topology (``.prmtop``) and coordinates (``.inpcrd``).

Step 5: Run a Simulation
------------------------

Run a short molecular dynamics simulation using OpenMM:

.. code-block:: python

   from simtk import openmm, unit
   from simtk.openmm import app

   pdb = app.PDBFile("pha_polymer.pdb")
   forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")

   system = forcefield.createSystem(pdb.topology,
                                    nonbondedMethod=app.PME,
                                    nonbondedCutoff=1.0*unit.nanometers,
                                    constraints=app.HBonds)

   integrator = openmm.LangevinIntegrator(
       300*unit.kelvin,   # Temperature
       1.0/unit.picoseconds, # Friction coefficient
       0.002*unit.picoseconds # Time step
   )

   simulation = app.Simulation(pdb.topology, system, integrator)
   simulation.context.setPositions(pdb.positions)

   simulation.minimizeEnergy()
   simulation.reporters.append(app.PDBReporter("pha_trajectory.pdb", 100))
   simulation.reporters.append(app.StateDataReporter("pha_log.txt", 100, step=True, potentialEnergy=True, temperature=True))

   simulation.step(1000)

This produces a trajectory file (``pha_trajectory.pdb``) and a log file (``pha_log.txt``).

Step 6: Analyze Results
-----------------------

Use Polymer Simulator’s analysis tools to process the trajectory:

.. code-block:: python

   from polymersimulator.analysis import Analyzer

   analyzer = Analyzer("pha_trajectory.pdb")
   analyzer.plot_rmsd("pha_rmsd.png")

This generates an RMSD plot of your PHA simulation.

Next Steps
----------

- Try different chain lengths (``n_units``).  
- Explore other simulation conditions (temperature, solvent models).  
- Check the :doc:`usage` section for more advanced workflows.

