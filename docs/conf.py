import os
import sys
sys.path.insert(0, os.path.abspath('../'))  # Adjusting to root project directory
sys.path.insert(0, os.path.abspath('../modules/'))  # Correcting to the directory of the modules

# Project information
project = 'Polymer Simulator'
copyright = '2025, Daniel J. York'
author = 'Daniel J. York'
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
html_static_path = ['_build/html']

master_doc = 'index'

# Autodoc default options (you can modify this as needed)
autodoc_default_options = {
    'members': True,            # Include members (functions, classes, etc.)
    'undoc-members': True,      # Include members without docstrings
    'show-inheritance': True,   # Show inheritance in class documentation
}
