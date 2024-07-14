# # Configuration file for the Sphinx documentation builder.
# #
# # For the full list of built-in configuration values, see the documentation:
# # https://www.sphinx-doc.org/en/master/usage/configuration.html

# import os
# import sys
# sys.path.insert(0, os.path.abspath('../../source'))

# # -- Project information -----------------------------------------------------
# # https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# project = 'SafeInCave'
# copyright = '2024, Herminio T. Honorio'
# author = 'Herminio T. Honorio'
# release = '0.1'

# # -- General configuration ---------------------------------------------------
# # https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# extensions = [
#     'sphinx.ext.autodoc',
#     'sphinx.ext.napoleon',
# ]

# templates_path = ['_templates']
# exclude_patterns = []



# # -- Options for HTML output -------------------------------------------------
# # https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
# html_static_path = ['_static']



# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Include source files ----------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../../safeincave'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SafeInCave'
copyright = '2024, Herminio T. Honorio'
author = 'Herminio T. Honorio'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    # 'sphinx.ext.viewcode',
    # 'sphinxcontrib.pseudocode',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_pdj_theme
html_theme_path = [sphinx_pdj_theme.get_html_theme_path()]

# html_theme = 'pydata_sphinx_theme'
# html_theme = 'piccolo_theme'
# html_theme = 'sphinx_pdj_theme'
html_theme = 'furo'
html_static_path = ['_static']

latex_elements = {
    'preamble': r'''
        \usepackage{tabularx}
        \usepackage{longtable}
        \usepackage{array}
        \usepackage{graphicx}
        \usepackage{fancyhdr}
        \fancypagestyle{plain}{
            \fancyhf{}
            \renewcommand{\headrulewidth}{0pt}
        }
        \pagestyle{plain}
    ''',
    'maketitle': r'''
        \begin{titlepage}
            \centering
            \vspace*{2cm}
            \Huge \textbf{Documentation} \\
            \vspace{1cm}
            \LARGE \textbf{Version 0.1.0} \\
            \vfill
            \includegraphics[width=0.3\textwidth]{logo_4.png}
            \vfill
            \Large \textbf{Herminio T. Honorio} \\
            \large \textbf{July 2024}
        \end{titlepage}
        \begingroup
          \pagestyle{empty}
          \cleardoublepage
        \endgroup
    ''',
}
