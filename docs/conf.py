import os
import sys
sys.path.insert(0, os.path.abspath('../'))  # Adjusting to root project directory
sys.path.insert(0, os.path.abspath('../modules/'))  # Correcting to the directory of the modules

from unittest import mock

# Mock openbabel to avoid build issues
sys.modules['openbabel'] = mock.MagicMock()
sys.modules['pybel'] = mock.MagicMock()


#from sw_openmm import *
#from sw_build_systems import *
#from sw_analysis import *
#from sw_basic_functions import *
#from sw_complex_fluid_models import *
#from sw_custom_decorators import *
#from sw_directories import *
#from sw_file_formatter import *
#from sw_depreceated import *
#from sw_orca import *

# Project information
project = 'Polymer Simulator'
copyright = '2025, Daniel J. York'
author = 'Daniel J. York, Isaac Vidal-Daza'
release = '0.1.0'

# Sphinx extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',  # Optional for Google/NumPy docstring support
]

# Paths and templates
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# HTML output
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    "collapse_navigation": False,  # Ensures the sidebar is expandable
    "sticky_navigation": True,     # Keeps the sidebar visible while scrolling
    "navigation_depth": 5,         # Allows deeper levels of nesting
}

html_static_path = ['_static']

master_doc = 'index'

# Autodoc default options (you can modify this as needed)
autodoc_default_options = {
    'members': True,            # Include members (functions, classes, etc.)
    'undoc-members': True,      # Include members without docstrings
    'show-inheritance': True,   # Show inheritance in class documentation
}
