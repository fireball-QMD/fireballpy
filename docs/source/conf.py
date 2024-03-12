# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))

import fireballpy

project = 'FireballPy'
copyright = '2024, José Ortega Mateo, Linda Angela Zotti, Jesús Ignacio Mendieta Moreno, Daniel González Trabada, Jorge Vega Martín, Carlos Roldán Piñero'
author = 'José Ortega Mateo, Linda Angela Zotti, Jesús Ignacio Mendieta Moreno, Daniel González Trabada, Jorge Vega Martín, Carlos Roldán Piñero'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'numpydoc',
    'sphinx_design',
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = '.rst'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static', os.path.join('..', '..', 'jupyter')]
