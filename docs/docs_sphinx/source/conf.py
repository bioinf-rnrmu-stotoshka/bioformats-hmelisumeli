import os
import sys

# путь от conf.py (docs_sphinx/source) до кода/code_scripts
sys.path.insert(0, os.path.abspath('../../code/code_scripts'))


html_theme = 'sphinx_rtd_theme'
# -- Project information -----------------------------------------------------
project = 'Infa_project'
copyright = '2025, Anastas'
author = 'Anastas'
release = '0.1'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',       # автогенерация из docstrings
    'sphinx.ext.napoleon',      # Google/NumPy стиль docstrings
    'sphinx.ext.viewcode',      # ссылки на исходники
]

templates_path = ['_templates']
exclude_patterns = []

language = 'ru'

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']



