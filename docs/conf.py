# -*- coding: utf-8 -*-
#
# hl7v2GenomicsExtractor documentation build configuration file, created by
# Rohan Gupta on Fri May 21 16:21:00 2021.
#
# This file is execfile()d with the current directory set to its containing
# dir.
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sphinx_rtd_theme
import sys
import os

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('..'))

# -- General configuration -----------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.intersphinx', 'sphinx.ext.autodoc']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
# source_encoding = 'utf-8'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'hl7v2GenomicsExtractor'
copyright = '2021, info@elimu.io'
author = ''


# The version info for the project you're documenting, acts as replacement for
# The short X.Y version.

# for source files.
exclude_trees = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output ---------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.

html_theme = "sphinx_rtd_theme"

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'collapse_navigation': False,

}

# mock depenpencies
autodoc_mock_imports = [
    "pyranges",
    "vcf",
    "hl7apy",
    "pytz",
    "pandas",
    "numpy"]

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []
html_static_path = ['_static']

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = ''

# Output file base name for HTML help builder.
htmlhelp_basename = 'hl7v2GenomicsExtractordoc'


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'http://docs.python.org/': None}
language = 'en'
# Configuration file for the Sphinx documentation builder.
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
