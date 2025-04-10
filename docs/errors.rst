Errors
======

Here some common errors that you may encounter are reported. Bear in mind that this seciton of the documentation will be updated over time
and some errors may not be logged yet and some may not have answers.

sw_openmm
=========

Periodic box size less than twice the nonbonded cutoff
------------------------------------------------------

This error relates the running of simulations using the **sw_openmm** module and this error will
throw an exception that will look similar to the error shown below. It usually occurs during an NPT simulation.

.. code-block:: bash

	openmm.OpenMMException: The periodic box size has decreased to less than twice the nonbonded cutoff.

This is a common and an expected error and is explained in the figure below.*(Note: of course a simulation takes place in 3D space, however this diagram
uses a 2D representation to showcase the issue)"

.. image:: images/PBC_error.png

The purple line indicates the nonbonded cutoff distance (this is set at 10 angstroms by default). This is the maximum distance where an atom can "see"
another atom in the simulation. During energy minimizations and NVT simulations the box will not change **but** during NPT simulations the box size
can grow//shrink. In cases where it shrinks **and** it shrinks to less than twice the nonbonded cutoff distance (i.e. less than 20 angstroms in any of 
the 3 directions {x,y,z}) the purple line will be outside of the box and an atom can potentially "see" itself in it periodic image. This only occurs 
if your system is small and equilibration takes your system to a small box or if you try and run an NPT simulation of a single molecule (standard box sizes
given to an individual molecule in this module are 20 angstroms in all directions).

Non parallel periodic box vectors
---------------------------------


This error relates the running of simulations using the **sw_openmm** module and this error will
throw an exception that will look similar to the error shown below.

.. code-block:: bash

	openmm.OpenMMException: First periodic box vector must be parallel to x.

This error implies that the periodic box vectors are not parallel and have deviated from the desired cubic shape. This error has not been seen before
and is being investigated.