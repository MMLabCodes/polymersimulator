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

**Antechamber** is part of AmberTools and carries out the parameterization of the molecules//polymers.

..literal-block::
	
	Antechamber

If antechamber is available, you will see something similar to the following:

.. image:: images/antechamber.PNG

**Tleap** is an AmberTools programme that allows for the building of systems and generation of topology and coordinate files for 

..literal-block::
	
	Tleap

If tleap is available, you will see something similar to the following:

.. image:: images/tleap.PNG

To exit tleap hold ctrl+c

A test for python and some associated packages is also important.

.. literal-block::
	
	python3

If python is available, you will see something similar to the following:

.. image:: images/python.PNG

*Note: this will open the python interpreter and you can enter any python commands after '>>>'*

With the python interpreter still open, check if the openmm python packages can be imported:

.. literal-block::
	
	from simtk.openmm import app

You may see the following warning:

.. image:: images/openmm_warning.PNG

This warning can be ignored, 'import openmm' is a better practice than the command above but importing from simtk will still load the openmm package.

If openmm is not installed properly you will see this:

.. image:: images/openmm_error.PNG

If you see this, exit the python interpreter with 'ctr+d' and execute the following in the command line:

.. literal-block::
	
	conda update -c conda-forge openmm

Now open the python interpreter again and try and import openmm.

Cloning the repository - Normal method
--------------------------------------

Now the repository needs to be cloned to give access to pyton modules within this pacakge. Go to the repository page in github https://github.com/MMLabCodes/polymersimulator and find the blue button labelled '<> code'.
Click here and select 'HTTPS' and copy the link. Now return to command line (ubuntu in windows), ensure you are in your home directory and execute the folowing:

..literal-block::
	
	git clone copied link

This will clone the repository into your pc and you will be able to access all the required files. 
Don't forget you can naviagate through the file explorer to view these files (if you are using windows, look for the linux folder with the penguin).

If this method doesn't work, see the alternative method below.


Cloning the repository - Alternative method
-------------------------------------------

First you will need to obtain a personal access token from github, once you have logged into github, click on your profile in the top right and navigate to (settings --> developer settings --> personal access tokens --> Tokens (classic)). 
Here, click on "generate new token --> generate new token (classic)" and enter a note "clone repo" and in the tick boxes, select "repo". 
Now scroll to the bottom and "generate token". This will give you a token you will need for the next step.

Now you can navigate to your home directory and execute the following commands:

..literal-block::
	
	git clone https://USERNAME:YOUR_TOKEN@github.com/MMLabCodes/polymersimulator.git
	cd polymersimulator

The final 'cd' command navigate to the directory containing the notebooks and scripts required for the tutorials.

Jupyter Notebook tutorials
==========================

There are a series of jupyter notebooks that contain tutorials and are aptly labelled. A section can be found explaining the contents of each notebook.
To launch jupyter notebooks execute the following in the command line (ubuntu in windows):

.. literal-block::
	
	jupyter notebook

This will launch a local jupyter notebook server and a series of URLs will be returned. Copy the **first link** containing ('localhost:8888') and copy//paste into a browser.
From there, refer to the section in this documentation about different tutorials.

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







