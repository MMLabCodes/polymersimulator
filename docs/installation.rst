Installation
============

This python package is designed to work within the unix shell. If you are working with **windows**, install miniconda on **ubuntu** https://ubuntu.com/ 

To get started you will need **anaconda** or **miniconda**:   

Anaconda: https://docs.anaconda.com/anaconda/install/   
Miniconda: https://docs.anaconda.com/miniconda/install/ (recommended for use with ubuntu)   

**Ubuntu** is a linux system emulator which is required to run some of the programmes utilised by this python package - please install miniconda using the ubuntu command line.

Creating an environment
=======================

To create an environment, open the command prompt (ubuntu if you are using windows) and execute the following:

..literal-block::

	conda create --name AmberTools23
	conda activate AmberTools23

*Note: You do not have to call your environment **AmberTools23** as long as the name makes sense to you*

Now the packages and different dependancies can be installed. These are all located in the **environment.yml** file in the **docs** folder of the github repository.
This file can be used to install of the required pacakges.

..literal-block::
	
	conda env update --file docs/environment.yml