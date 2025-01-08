Installation
============

This python package is designed to work within the unix shell. If you are working with **windows**, install miniconda on **ubuntu** https://ubuntu.com/ 

To get started you will need **anaconda** or **miniconda**:   

Anaconda: https://docs.anaconda.com/anaconda/install/   
Miniconda: https://docs.anaconda.com/miniconda/install/ (recommended for use with ubuntu)   

**Ubuntu** is a linux system emulator which is required to run some of the programmes utilised by this python package - please install miniconda using the ubuntu command line.

Creating an environment
-----------------------

To create an environment, open the command prompt (ubuntu if you are using windows) and execute the following:

..literal-block::

	conda create --name AmberTools23
	conda activate AmberTools23

*Note: You do not have to call your environment **AmberTools23** as long as the name makes sense to you*

Now the packages and different dependancies can be installed. These are all located in the **environment.yml** file in the **docs** folder of the github repository.
This file can be used to install of the required pacakges.

..literal-block::
	
	conda env update --file docs/environment.yml

Testing
-------

To ensure some key packages are installed it is important to run some tests. Execute the commands below and compare with the images.

Antechamber is part of AmberTools and carries out the parameterization of the molecules//polymers.

..literal-block::
	
	Antechamber

If antechamber is available, you will see something similar to the following:

.. image:: images/antechamber.png

Tleap is an AmberTools programme that allows for the building of systems and generation of topology and coordinate files for 

..literal-block::
	
	Tleap

If tleap is available, you will see something similar to the following:

.. image:: images/tleap.png

To exit tleap hold ctrl+c

A test for python and some associated packages is also important.

..literal-block::
	
	python3

If python is available, you will see something similar to the following:

.. image:: images/python.png

*Note: this will open the python interpreter and you can enter any python commands after '>>>'*

With the python interpreter still open, check if the openmm python packages can be imported:

..literal-block::
	
	from simtk.openmm import app

You may see the following warning:

.. image:: images/openmm_warning.png

This warning can be ignored, 'import openmm' is a better practice than the command above but importing from simtk will still load the openmm package.

If openmm is not installed properly you will see this:

.. image:: images/openmm_error.png

If you see this, exit the python interpreter with 'ctr+d' and execute the following in the command line:

..literal-block::
	
	conda update -c conda-forge openmm

Now open the python interpreter again and try and import openmm.






