import os
import sys
sys.path.insert(0, os.path.abspath('../'))  # Adjusting to root project directory

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
html_static_path = ['_static']

master_doc = 'index'

autodic_default_options = {'members':True, 'undoc-members':True, 'show-inheritance':True}
