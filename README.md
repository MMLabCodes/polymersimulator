# Polymer Simulator

This repository contains content that will allow building arrays of polymer chains, parametrise them with **Amber** force field, and perform simulations using **openmm**

## 1. Setting up the environment and cloning repository

(Important: If you are working in windows, please follow the steps in section section 4 <a name="section-4"></a> (Working with windows) to first set up a linux system in your computer before following the steps to set up the environment. If you are working in Unix (linux or macOS), you will not require any prerequisite steps before setting up the environment given below)

To run this code, a Python environment containing **RDkit**, **AmberTools** and **openmm** is required. This environment is set up by running the lines described in the following steps 1-4 in your terminal (command line). (In windows, open ubuntu and enter these lines into the terminal)

1. Install miniconda

   Miniconda is a package and environment manager for Python. The following commands install miniconda in your home directory and intialise it for setting up environments.
   
```
cd 
mkdir -p ~/miniconda3 
cd miniconda3/ 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh 
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
~/miniconda3/bin/conda init bash 
~/miniconda3/bin/conda init zsh
```
2. Create environment

Creating an environment is useful as it creates a separate container for all packages required in the project. 
The following commands create an environment called "AmberTools23" and activate it.
   ```
   conda create --name AmberTools23 
   conda activate AmberTools23 
   ```
3. Install **openmm**, **RDkit** and **AmberTools**

   Now the environment has been activated the following commands install the 3 packages required.
   
   ```
   conda install -c conda-forge ambertools=23
   conda install -c conda-forge openmm
   conda install -c conda-forge rdkit
   conda install -c conda-forge openbabel
   ```
   
3.1 Install **Jupyter notebooks**

   We need to install Jupyter notebooks in this environment to run AmberTools from a Python notebook (in addition to running it from the terminal)
    
   ```
   sudo apt install python3-pip python3-dev
   pip install jupyter
   ```
    
3.2 Test the Jupyter notebook

   In the terminal (or Ubuntu terminal if using a Windows machine), we need to enter the following:
    
   ```
   jupyter notebook
   ```
   This will start a remote Jupyter notebook server within the environment we have just set up.
    
   You will see the following prompt after entering 'Jupyter notebook':

<img width="510" alt="jupyter_tut" src="https://github.com/DanielYyork/polymer_simulator/assets/93723782/9718b875-aeb0-421b-a134-87e945d9b585">

   Now, you can select the first URL (the one containing 'localhost:8888') and copy and paste it into a browser, this will launch Jupyter notebook (fingers crossed!)
    
   From there you can navigate to the Jupyter notebook folder and launch any additional notebooks from there.
    
   For now, we will close the notebook and ensure our other packages are working properly.
    
   To close the notebook, return to Ubuntu: hold "ctr" + "c" at the same time, you will be asked if you want to close Jupyter notebook - yes!
    
4. Ensure the required packages are available

   Before running any code, it is recommended to check the availability of different packages.

   **AmberTools** is a collection of different tools. Among them, **antechamber** and **tleap** are used extensively. To check these tools are available, enter the commands below into your terminal.
   ```
   antechamber
   ```
   If antechamber is available, you will see this in your terminal:
   ```
   Welcome to antechamber 22.0: molecular input file processor.
   Usage: antechamber -i     input file name
                   -fi    input file format
                   -o     output file name
                   -fo    output file format
                   -c     charge method
                   -cf    charge file name
                   -nc    net molecular charge (int)
                   -a     additional file name
                   -fa    additional file format
                   -ao    additional file operation
                   ... more operations ..
   ```

   ```
   tleap
   ```
   If tleap is available, you will see this in your terminal:
   ```
   Welcome to LEaP!
   (no leaprc in search path)
   >
   ```
   Note: The 'tleap' command opens an interactive version of the tool where you can enter tleap commands. To exit this press ctr+c.

   Checking the availability of **openmm** is slightly different as it is a Python package - not a standalone program. Open the Python interpreter as follows:
   ```
   python3
   ```
   If python3 is available, you will see this in your terminal:
   ```
   Python 3.12.1 | packaged by conda-forge | (main, Dec 23 2023, 08:03:24) [GCC 12.3.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>>
   ```
   Like 'tleap' the 'python3' command opens an interactive version of Python and Python code can be entered after '>>>'. To check openmm is imported enter the following into the Python interpreter:
   ```
   >>> from simtk.openmm import app
   ```
   You may see the following warning:
   ```
   Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
   ```
   This warning can be ignored, "import openmm" is better but importing from simtk will still load the openmm package.
   Note: This is an interactive Python interpreter and pressing ctr+d will exit this - ctr+c acts as a keyboard interrupt in the python interpreter and will interrupt any running code but will not exit.

   If openmm is not installed properly you will see this:
   ```
   >>> from simtk.openmm import app
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   ModuleNotFoundError: No module named 'simtk'
   >>>
   ```
   In this case return to step 3 and try to install openmm again. An update to the openmm package may also be the issue; it can be updated with the following line:
   ```
   conda update -c conda-forge openmm
   ```
   Now, try an import openmm to the python interpreter again - it should work!

5. Cloning the repository

5.1 Normal git clone method
    
   At the top of the GitHub page there will a blue button labelled '<> code'. Click here and select 'HTTPS' and copy the link. No return to Ubuntu and enter:
    
   ```
   git clone copied_link
   ``` 
   This will clone the repository into Ubuntu and you will be able to access all the required files.
   Don't forget you can navigate through the file explorer to view these files (see section 4 for where linux files are located).
    
   If this method does not work, see the alternative method below.

5.2 Alternative git clone method

   To download these Python scripts and Jupyter notebooks it is necessary to clone the BCSW repository. **You will need a GitHub account**.
   This will give you access to all the files in your own computer. The commands below should be executed in a terminal,
   this will create a new directory in you home directory.
    
   First you will need to obtain a personal access token from GitHub, once you have logged into GitHub, 
   click on your profile in the top right and navigate to (settings --> developer settings --> personal access tokens --> Tokens (classic)). 
   Here, click on "generate new token --> generate new token (classic)" and enter a note "clone repo" and in the tick boxes, select "repo". 
   Now scroll to the bottom and "generate token".
   This will give you a token you will need for the next step.
   ```
   cd 
   git clone https://USERNAME:YOUR_TOKEN@github.com/MMLabCodes/BCSW.git
   cd polymer_simulator
   ```
   The final command navigates to the directory containing the notebooks and scripts required for the tutorial.

## 2. Jupyter notebooks

There are a series of Jupyter notebooks that contain examples and detailed explanations of how to carry out several tasks. In principle, 
these Jupyter notebooks are guides and explanations.
You can create your own notebook to combine the building of a system and analysis into one using examples from the prerequisite notebooks.
Here, short descriptions of each notebook are detailed. They are in order in which they should be viewed/executed.

We can launch Jupyter notebooks from Ubuntu with:

   ```
   jupyter notebook
   ```
No we can select the first URL (the one containing 'localhost:8888') and copy and paste it into a browser, this will launch Jupyter notebook (fingers crossed!)
 
From there we can navigate to the Jupyter notebook folder and launch notebooks from there.

### 2.1 Jupyter notebook tutorials

#### 2.1.1 Filepath manager

This tutorial provides examples of how to use the filepath manager. There will be many cases where nothing is returned as no
files are currently generated. However, it provides a useful framework for working with files and systems easily. Many of
the other methods described in future notebooks use this filepath manager so it is crucial to understand its functionality.

Do the tutorial: **Tutorial_1_filepath_manager.ipynb**

#### 2.1.2 Parameterizing small molecules and polymers

This notebook is split into 2 sections; <br>
    1. Parameterizing small molecules <br>
    2. Parameterizing polymers <br>
    
Parameterizing molecules is carried out using the code in **modules/sw_build_systems** and utilises amber tools to generate 
parameters using GAFF (generalized amber forcefield).

Do the tutorial: **Tutorial_2_Parameterizing_Small_Molecules_and_Polymers.ipynb**

*Note: this code for generating parameters is project agnostic, and for specific cases - easier tools may alreadt exist*

#### 2.1.3 Solvating small molecules and polymers

This notebook is split into 2 sections; <br>
    1. Solvating small molecules <br>
    2. Solvating polymers <br>
    
Amber topologies and parameters will be generated for systems of solvated small molecules and polymers to be used for molecular dynamics simulations. <br>

Do the tutorial: **Tutorial_3_Solvating_Small_Molecules_and_Polymers.ipynb**

#### 2.1.4 Building systems with polymers

This notebook is split into 2 sections; <br>
    1. Building 3x3 arrays of polymers <br>
    2. Building 2x10 arrays of polymers - *array is generated, but solvation of the array in water does not* <br>
    
Amber topologies and parameters will be generated for systesms of the above systesm to be used for molecular dynamics simulations. <br>

#### 2.1.5 Running simulations

This notebook will explain how to run simulations on systems where amber topologies and coordinates have been generated for. The example in this tutorial is for the 3x3 array of 3HB_10-mers generated in **tutorial_4**, however, this will work with any system generated in either **tutorial_3** or **tutorial_4**.


#### 2.2.2.2 Module guides

There are a series of guides that give a detailed over of what each module does - some of this information may also be included within the tutorials.

1. **openmm_simulation_guide.ipynb** <br>

    This notebook contains detailed explanations and examples of how to use the functions found in the **sw_openmm.py** module. This guide is essential to using this python module and should be read and referred to when running openmm simulations. This module functions as wrapper for openmm scripts and no experience with openmm is required to use this module. <br>

2. **parameterization_guide.ipynb** <br>

   This notebook includes details on parameterizing molecules for MD simulations using tleap and antechamber. All of methods shown here are included in the **sw_build_systems.py** module - it is recommend to work through the tutroial notebooks first and then refer to this notebook if you would like more information. <br>

3. **building_systems_guide.ipynb** <br>

    Notebook and associated description coming soon...

3. **filepath_manager_guide.ipynb** <br>

    Notebook and associated description coming soon...

### 3 Python modules - easy handling 

The python modules created here aim to provide easier handling of a range of different aspects of this project and are detailed in more depth below. <br>

1. **sw_directories**

    This module defines the PolymerSimulatorDirs class, which is used for initializing
    and organizing directories for polymer simulation setup. It creates specific
    directories such as 'python_scripts', 'pdb_files', 'molecules', 'systems', and manages
    files like 'residue_codes.csv'. More information can be found in the **filepath_manager_guide.ipynb** notebook.

    Attributes:
        main_dir (str): The main directory for polymer simulation setup.

    Example:
    ```
    >>> from PolymerSimulatorDirs import PolymerSimulatorDirs 
    >>> polymer_dirs = PolymerSimulatorDirs('/path/to/main/dir/')
    >>> print(polymer_dirs.pdb_file_dir)
    '/path/to/main/dir/pdb_files'
    ```

    Note:
        The main_dir must be provided as a valid path and must already exist.

2. **sw_openmm.py**

   This module defines the following classes:
   - BuildSimulation
   - AmberSimulation
   - ANISimulation <br>

   Functions for carrying MD simulations with openmm are contrained here; such as - annealing, equilibration and production runs. The jupyer notebook **openmm_simulation_guide.ipynb** contains examples and decriptions on how to use this module.

3. **sw_build_systems.py**

   Description coming soon...

4. **sw_analysis.py**

   Decription coming soon...


## 4. Working on windows system <a name="section-4"></a>

The Amber package is made for linux and it is best to install a linux distribution for windows following the instructions below. This will allow you run all of the code here in these notebooks and python scripts in a linux environment. One thing to keep in mind is that any simulations run in your local pc are simply tests before moving to Super Computing Wales (SCW).

(Important: if you are working on a linux or macos systemm - please skip to "setting up environment")

1. Enable windows subsystem for linux.

   In windows powershell (start menu --> windows powershell --> run as administrator) entering the following will enable
   windows subsytem for linux.
   ```
   dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
   ```
2. Install linux distribution - Ubuntu

   Ubuntu is a linux distribution and is available to download from the microsoft store:
   https://apps.microsoft.com/search?query=ubuntu&hl=en-us&gl=US

3. Launch ubuntu

   (start menu --> unbuntu) When you first launch ubuntu you will be asked to enter a username as password. This can whatever you desire but will what is used to log into to ubuntu in the future.

   You should now see something similar to this:
   ![image](https://github.com/MMLabCodes/BCSW/assets/93723782/fa03d5cf-cddc-4f87-97f8-85589e94b44a)

4. Upgrade  ubuntu packages

   ```
   sudo apt update
   sudo apt upgrade
   ```

   Note: you will be able to see linux in your file explorer.
   ![image](https://github.com/MMLabCodes/BCSW/assets/93723782/1c11cfe9-0182-4755-a0c9-06395fbd9767)

Now you can follow the instructions detailed above and run all python scripts. <br>

5. Installing and running Jupyter notebooks with Ubuntu.

Launching Jupyter notebooks is also possible from command line using linux and macOS with the created environment - a solution for doing this using Ubuntu is coming.



   


