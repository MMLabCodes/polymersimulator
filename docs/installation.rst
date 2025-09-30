Installation
============

The **Polymer Simulator** package is designed to work in a Unix-like shell.  
If you are on **Windows**, we recommend installing **Miniconda** inside **Ubuntu** (via WSL) to ensure compatibility:

- Ubuntu: https://ubuntu.com/  
- Anaconda: https://docs.anaconda.com/anaconda/install/  
- Miniconda (recommended): https://docs.anaconda.com/miniconda/install/

.. note::
   Ubuntu (via WSL on Windows) is required because some external programs used by this package are only available in Linux environments.

Cloning the repository
----------------------

Clone the GitHub repository to access the Python modules.

**Method 1: Standard Clone**

1. Go to the repo: https://github.com/MMLabCodes/polymersimulator  
2. Click the blue **<> Code** button, select **HTTPS**, and copy the link.  
3. From your terminal (Ubuntu/WSL), run:

.. code-block:: bash

   git clone <copied-link>

This will clone the repository into your home directory.

**Method 2: Using a Personal Access Token**

If the standard method fails, you can use a GitHub personal access token:

1. In GitHub, go to:  
   ``Settings → Developer settings → Personal access tokens → Tokens (classic)``  
2. Generate a new token with the ``repo`` permission.  
3. Use the token to clone:

.. code-block:: bash

   git clone https://USERNAME:YOUR_TOKEN@github.com/MMLabCodes/polymersimulator.git
   cd polymersimulator

The final ``cd`` command places you inside the repo where tutorials and scripts are stored.

Creating a Conda Environment
----------------------------

First, create and activate a new Conda environment:

.. code-block:: bash

   conda create --name polymer_simulator
   conda activate polymer_simulator

.. tip::
   The environment name does not have to be named ``polymer_simulator``. Choose any name that makes sense to you.

Next, install all required dependencies using the provided ``environment.yml`` file:

.. code-block:: bash

   conda env update --file docs/environment.yml

Testing the Installation
------------------------

To confirm everything is installed correctly, run the following checks.

**1. AmberTools (Antechamber)**

.. code-block:: bash

   antechamber

If successful, you should see output similar to:

.. image:: images/antechamber.PNG

**2. AmberTools (Tleap)**

.. code-block:: bash

   tleap

If successful, you should see output similar to:

.. image:: images/tleap.PNG

To exit ``tleap``, press ``Ctrl+C``.

**3. Python and OpenMM**

Check Python is installed:

.. code-block:: bash

   python3

You should see something like:

.. image:: images/python.PNG

Now test importing OpenMM:

.. code-block:: python

   from simtk.openmm import app

.. warning::
   You may see a warning about ``simtk``. This can be ignored, but using ``import openmm`` is the preferred practice.

If OpenMM is not installed correctly, you may see an error:

.. image:: images/openmm_error.PNG

To fix this, update OpenMM:

.. code-block:: bash

   conda update -c conda-forge openmm

Then retry the import.


Launching Jupyter Notebook
--------------------------

The repository includes Jupyter Notebook tutorials. You can launch them from your linux command prompt:

.. code-block:: bash

   jupyter notebook

This will start a local server and display several URLs. Copy the **first link** (containing ``localhost:8888``) into your browser.

From there, explore the tutorials or follow the detailed guides in this documentation.

