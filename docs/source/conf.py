# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import re
import fireballpy

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# General information
with open('../../pyproject.toml', 'rb') as fp:
    pyproject = tomllib.load(fp)
author = ', '.join([x['name'] for x in pyproject['project']['authors']])
project = 'FireballPy'
copyright = f'2024, {author}'
version = re.sub(r'\.dev.*$', r'.dev', fireballpy.__version__)
release = version

# Sphinx
extensions = [
    # 'sphinx_rtd_theme',
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'numpydoc',
    'sphinx_design',
    'sphinx_favicon',
    "sphinx_copybutton",
    'sphinxcontrib.bibtex', 
]

bibtex_bibfiles = ['references.bib']

master_doc = 'index'
source_suffix = '.rst'

today_fmt = '%B %d, %Y'
templates_path = ['_templates']
exclude_patterns = []

# HTML
# html_theme = 'sphinx_rtd_theme'
html_theme = 'pydata_sphinx_theme'
html_logo = '_static/fireball.svg'
html_favicon = '_static/fireball.ico'
html_sidebars = {
    "index": ["search-button-field"],
    "**": ["search-button-field", "sidebar-nav-bs"]
}
html_theme_options = {
    "external_links": [
        {
            "url": "https://fireball-qmd.github.io",
            "name": "Fireball",
        },
        {
            "url": "https://wiki.fysik.dtu.dk/ase/",
            "name": "ASE",
        },
    ],
    "header_links_before_dropdown": 4,
    "logo": {
        "text": "FireballPy",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/fireball-QMD/fireballpy",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/fireballpy/",
            "icon": "fa-custom fa-pypi",
        },
    ],
    "navbar_start": ["navbar-logo"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_persistent": [],
    "secondary_sidebar_items": ["page-toc"],
    "footer_start": ["copyright"],
    "footer_center": ["sphinx-version"],
}

html_title = f"{project} v{version} Manual"
html_static_path = ['_static', os.path.join('..', '..', 'examples')]
html_last_updated_fmt = '%b %d, %Y'
html_additional_pages = {}
html_use_modindex = True
html_domain_indices = False
html_copy_source = False
html_file_suffix = '.html'
htmlhelp_basename = 'fireballpy'
html_js_files = ["custom-icon.js"]
html_css_files = ["fireballpy.css"]

# Copy-button
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.{3,}: | {5,8}: "
copybutton_prompt_is_regexp = True

# Autosummary
autosummary_generate = True

# Autodoc
autodoc_default_options = {
    'inherited-members': None,
}
autodoc_typehints = 'none'
