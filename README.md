# Polymer simulator

This directory contains content that will allow for simulations of polymers using **openmm** with **Amber**.

## 1. Setting up the environment and cloning repository

(Important: If you are working in windows, please follow the steps in section section 4 <a name="section-4"></a> (Working with windows) to first set up a linux distribution in your computer before following the steps to set up the environment. If you are working in linux or macos, you will not require any prerequisite steps to setting up the environment given below)

To run the code in this folder, an environment containing **RDkit**, **AmberTools** and **openmm** is required and is set up by running the lines steps 1-4 in your command line. (In windows, open ubuntu and enter these lines into that terminal)

1. Install miniconda

   Miniconda is a package and environment manager for python. The following commands install miniconda in your home directory and intialise for setting up environments.
   
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

Creating an environment is useful as it creates a seperated container for all of the packages required in this project. 
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
   ```
4. Ensure packages are available

   Before running any code, it is recommended to check the availiability of different packages.

   **AmberTools** is a collection of different programmes and 2 programmes used extensivley are **antechamber** and **tleap**. To check these are available, enter the commands below into your command line.
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
   Note: The 'tleap' command opens an interactive version of the programme where you can enter tleap commands. To exit this press ctr+c silmultaneously.

   Checking if **openmm** is available is slightly different as it a python package - not a standalone programme. open the python interpreter as follows:
   ```
   python3
   ```
   If python3 is available, you will see this in your terminal:
   ```
   Python 3.12.1 | packaged by conda-forge | (main, Dec 23 2023, 08:03:24) [GCC 12.3.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>>
   ```
   Like 'tleap' the 'python3' command opens an interactive version of python and python code can be entered after '>>>'. To check openmm is imported enter the following into the python interpreter:
   ```
   >>> from simtk.openmm import app
   ```
   You may see the following warning:
   ```
   Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
   ```
   This warning can be ignored, "import openmm" is better but importing from simtk will still load the openmm package.
   Note: This is an interactive python interpreter and pressing ctr+d silmulataneously will exit this - ctr+c acts as a keyboard interrupt in the python interpretor and will interupt any running code but will not exit.

   If openmm is not installed properly you will see this:
   ```
   >>> from simtk.openmm import app
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   ModuleNotFoundError: No module named 'simtk'
   >>>
   ```
   In this case return to step 3 and try to install openmm again. An update to the openmm package may also be the issue and it can be updated with the following line:
   ```
   conda update -c conda-forge openmm
   ```
   Now try an import openmm to the python interpreter again - it should work!

5. Cloning the repository

    To download these python scripts and jupyter notebooks it is necessary to clone the BCSW repository. **You will need a github account**.
    This will give you access to all of the files in your own computer. The commands below should be executed in command line (linux/macos) or
    ubuntu if you are using windows - this will create a new directory in you home directory.
    
    First you will need to obtain a personal access token from github, once you have logged into github, click on your profile in the top right and navigate to 
    (settings --> developer settings --> personal access tokens --> Tokens (classic)). Here, click on "generate new token --> generate new token (classic)" and enter a note "clone repo" and in the tick boxes, select "repo". Now scroll to the bottom and "generate token".
    This will give you a token you will need for the next step.
    ```
    cd 
    git clone https://USERNAME:YOUR_TOKEN@github.com/MMLabCodes/BCSW.git
    cd BCSW/examples/openmm/polymer_simulator
    ```
    The final command navigates to the directory containing the notebooks and scripts required for the tutorial.
    
6. Running Jupyter notebooks from ubuntu

    Issues encountered - it is easy enough to launch from anaconda and run independently to ubuntu, but there is a way... solution coming soon.

## 2. Jupyter notebooks

IMPORTANT: Notebooks currently can't be launched from ubuntu, but this will be fixed soon. It should allow for easy execution of both command line and python orientated notebooks.

There are a series of jupyter notebooks that contain examples and detailed explanations of how to carry out several tasks. In principal these jupyter notebooks are guides and explanations and there will be an associated python script.
Here, short descriptions of each notebook are detailed and there is an order in which to view/execute them.

1. **parametrization.ipynb**

    This notebook contains information on how to parametrize molecules from pdb files and uses the 3HB monomer as an example. 
    The inputs folder contains the 3HB.pdb file.

2. **python_parametrization.ipynb**

    Description coming soon - needs to be integrated with build_simulation functions (i.e. directory creation and generation of filepaths)

3. **openmm_simulations.ipynb**

    This notebook contains examples and explanations of steps to take in openmm to run simulations.
    This includes information about different simulations:
        1. Simple MD simulation (production run)
        2. Equilibration of a system followed by a production run
            (Equilibration is where a system reaches a stable state at the correct density)

## 3. Python scripts 

There are a series of python scripts associated with this project:
    - pdb file generation 
    - Molecular dynamics input (paramater/topology) generation
    - Molecular dynamics scripts

~~Descriptions coming soon...~~
Decsriptions currently being added...

### 3.1 Python modules - easy handling 

The python modules created here aim to provide easier handling of a range of different aspects of this project and are detailed in more depth below. <br>

1. **sw_directories**

    This module defines the PolymerSimulatorDirs class, which is used for initializing
    and organizing directories for polymer simulation setup. It creates specific
    directories such as 'python_scripts', 'pdb_files', 'molecules', 'systems', and manages
    files like 'residue_codes.csv'.

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

### 3.2 Python scripts - pdb file generation

1. **smiles_to_pdb.py**

    This python script generates a PDB file from a SMILES string and also generates a unique residue code
        for that PDB file by generating and cross-referencing with a database of all generated pdbs. This residue code
        is useful for tracking and analysis of molecules in MD simulations.

    This python file requires the following inputs: <br>
    ```
    Parameters: 
    - smiles (str): SMILES representation of the molecule. 
    - name (str): Name of the molecule. 
    - directory (str): Directory where the PDB file will be saved. 
    - residue_code_csv (str): Path to the CSV file containing existing residue codes. 
    ```
        
    This python script should be run from the "polymer_simulator" directory as follows: <br>
    ```
    python3 python_scripts/smiles_to_pdb.py C methane /pdb_files/molecules /pdb_files/residue_codes.csv 
    ```
    The general form of providing arguments for this script is: <br>
    ```
    python3 python_scripts/smiles_to_pdb.py <smiles> <molecule_name> <pdb_file_directory> <residue_code_csv>
    ```
    
    It is necessary to create the "pdb_files/molecules" directory intially, but the csv should be generated automatically if it does not already exist.

2. **csv_to_pdb.py**

Coming soon...

### 3.2 Python scripts - Molecular dynamics input (paramater/topology) generation 
    
1. **python_parametrizer.py**

    For this section of python scripts - the output from this script is the most important 
        and is required by each subsequent script. In short, this python script parametrizes
        the molecules and the subsequent scripts build systems by calling these parameters. 
        
    This python file requires the following inputs:
    ```
    Parameters:
    - pdb file
    ```
    
    This python script should be run the the "polymer_simulator" directory as follows :
    ```
    python3 python_scripts/python_parametrizer.py path_to_pdb_file.pdb 
    ```
    The general form of providing arguments for this script is: <br>
    ```
    python3 python_scripts/python_parametrizer.py <path_to_pdb_file.pdb> 
    ```
    The "path_to_pdb_file.pdb" path should be something similar: <br>
        "~/polymer_simulator/pdb_files/molecules/pdb_file.pdb"

2. **single_molecule_solvater.py**

    This python script generates the topology and paramaters for a system comprised of a single solvated
        monomer using *AMBER* for molecular dynamics simulations. 
        
    This python file requires the following inputs: <br>
    ```
    Parameters: 
    - molecule_name (str): Name of the molecule. 
    ```
    This python script should be run from the "polymer_simulator" directory as follows: <br>
    ```
    python3 python_scripts/single_molecule_solvater.py methane 
    ```  
    The general form of providing arguments for this script is: <br>
    ```   
    python3 python_scripts/single_molecule_solvater.py <molecule_name> 
    ```
    The AMBER files generated are for a single solvated molecule in a box with dimensions twice 
        that of its max pairwise distance.
        
    If the molecule being solvated has not been previously solvated with **python_parametrizer.py** and
        error will be seen asking you to first parametrize the molecue.

6. **2_2_array_generator.py**

Coming soon...

7. **3_3_array_generator.py**

Coming soon...

### 3.3 Python scripts - Molecular dynamics scripts 

ADD ABOUT GENERATED INPUTS BEING USED FOR THE SCRIPTS - ADD A PART ABOUIT WHAT INPUTS ARE GENERATED FROM EACH SCRIPT
MAYBE ADD PICS OF THEM AS WELL <br>

The prmtop and inpcrd files are generated by the scripts given in section "3.2 python scripts for molecular dynamics input"
and each python script in that section will generate valid inputs for this script. <br>
        
*Note: Other types of paramater and coordiante files will be supported, but for now only **".prmtop"** and **".rst7"** filetypes are supported.*

1. **simple_md_run.py**

    This python script carries out a simple molecular dynamics run - an NVT simulation.
    The default length of this simulation is 10,000 steps with a timestep fo 2 fs - these can be changed manually by editting the script.

    This python file requires the following inputs: <br>
    ```
    Parameters: 
    - prmtop file (str): Parameter file for a system (this file will have the ".prmtop" extension) 
    - inpcrd file (str): Coordinate file for a system (this file will have the ".rst7" extension) 
    ```
    This python script should be run from the "polymer_simulator" directory as follows: <br>
    ```
    python3 python_script/simple_md_run.py path_to_prmtopfile path_to_rst7file 
    ```
    The general form of providing arguments for this script is: <br>
    ```  
    python3 python_scripts/simple_md_run.py <parameter_file> <coordinate_file> 
    ```

      
2. **equilibration_production.py**

    This python script carries out a molecular dynamics simulation in 2 parts: <br>
        1. Equilibration for 1,000,000 steps with a timestep of 2fs (NPT ensemble) <br>
        2. Production run for 2,000,000 steps with a timesept of 2fs (NVT ensemble) <br>
        
    This python file requires the following inputs: <br>
    ```
    Parameters: 
    - prmtop file (str): Parameter file for a system (this file will have the ".prmtop" extension) 
    - inpcrd file (str): Coordinate file for a system (this file will have the ".rst7" extension) 
    ```
        
    This python script should be run from the "polymer_simulator" directory as follows: <br>
    ```
    python3 python_scripts/equilibration_production.py path_to_prmtopfile path_to_rst7file <br>
    ```   
    The general form of providing arguments for this script is: 
    ```  
    python3 python_scripts/equilibration_production.py <parameter_file> <coordinate_file> 
    ```   

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

5. Installing and running jupyter notebooks with ubuntu.

Launching jupyter notebooks is also possilbe from command line using linux and macos with the created environment - a solution for doing this using ubuntu is coming.



   


