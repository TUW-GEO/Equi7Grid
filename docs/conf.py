# noqa: INP001
"""Sphinx configuration for generating documentation.

This file only contains a selection of the most common options.
For a full list see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html

"""

import sys
from pathlib import Path

__location__ = Path(__file__).parent
sys.path.insert(0, (__location__ / "../src").resolve().as_posix())

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "myst_nb",
    "autodoc2",
    "sphinx.ext.napoleon",
]
myst_heading_anchors = 3
autodoc2_packages = [
    "../src/equi7grid",  # path to your package
]
autodoc2_render_plugin = "myst"
autodoc2_hidden_objects = ["dunder", "undoc", "private", "inherited"]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix of source filenames.
source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "equi7grid"
copyright = "2025, TU Wien"  # noqa: A001

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".venv"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"


# If this is True, todo emits a warning for each TODO entries. The default is False.
todo_emit_warnings = True

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_book_theme"

# Output file base name for HTML help builder.
htmlhelp_basename = "equi7grid-doc"

html_theme_options = {
    "switcher": {
        "json_url": "https://tuw-geo.github.io/equi7grid/switcher.json",
        "version_match": "latest",
    },
    "navbar_end": ["version-switcher"],
}
