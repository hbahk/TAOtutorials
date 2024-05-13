# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TAO Tutorials'
copyright = '2024, Hyeonguk Bahk'
author = 'Hyeonguk Bahk'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "myst_nb",
    "sphinx_togglebutton",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "sphinx_sitemap",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
    '.myst': 'myst-nb',
}

myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

jupyter_execute_notebooks = "off"

nb_render_image_options = {"align": "center"}
    
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

html_css_files = [
    'custom.css',
]

html_theme_options = {
    "repository_url": "https://github.com/hbahk/TAOTutorials",
    "use_repository_button": True,
    "logo": {
        "text": "TAO Tutorials",
        "alt_text": "TAO Tutorials - Home",
    }
}

html_baseurl = "https://hbahk.github.io/TAOTutorials/"

html_js_files = ["https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js"]