# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------


import os
import re
import subprocess
import sys
from pathlib import Path

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "ClassyQ"
copyright = "2022, Roberto Di Remigio Eikås and contributors"
author = "Roberto Di Remigio Eikås"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.autosectionlabel",
    "sphinxcontrib.bibtex",
    "breathe",
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Extentions configuration ------------------------------------------------

todo_include_todos = True

breathe_projects = {"ClassyQ": "_build/xml"}
breathe_default_project = "ClassyQ"
breathe_default_members = ("members", "protected-members", "private-members")

bibtex_bibfiles = ["bibliography.bib"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_theme_options = {
    "light_logo": "ClassyQ.png",
    "dark_logo": "ClassyQ.png",
}


major = "0"
minor = "0"
patch = "0"

version = f"{major}.{minor}.{patch}"


class SearchReplace(dict):
    """All-in-one multiple-string-substitution class."""

    def _make_regex(self):
        """Build re object based on the keys of the current dictionary."""
        return re.compile("|".join(map(re.escape, self.keys())))

    def __call__(self, match):
        """Handler invoked for each regex match."""
        return self[match.group(0)]

    def replace(self, text):
        """Translate text, returns the modified text."""
        return self._make_regex().sub(self, text)


def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""

    try:
        retcode = subprocess.call(f"cd {folder}; doxygen", shell=True)
        if retcode < 0:
            sys.stderr.write(f"doxygen terminated with signal {retcode}")
    except OSError as e:
        sys.stderr.write(f"doxygen execution failed: {e}")


def setup(app):

    # are we on the readthedocs servers?
    read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"

    # project_root_dir -- the root of the project
    # project_src_dir  -- source code location
    # project_docs_dir  -- .rst location
    project_docs_dir = Path.cwd()
    project_root_dir = project_docs_dir.parent
    project_src_dir = project_root_dir / "src"
    print(f"Project root directory {project_root_dir}")
    print(f"Project doc directory {project_docs_dir}")
    print(f"Project src directory {project_src_dir}")

    # configure Doxyfile.in
    rep = {
        "@PROJECT_VERSION_MAJOR@": major,
        "@PROJECT_VERSION_MINOR@": minor,
        "@PROJECT_VERSION_PATCH@": patch,
        "@PROJECT_SOURCE_DIR@": str(project_root_dir),
    }
    replacer = SearchReplace(rep)

    # read in docs/Doxyfile.in
    doxyfile_in = project_docs_dir / "Doxyfile.in"
    doxyfile = project_docs_dir / "Doxyfile"

    with doxyfile_in.open("r") as f:
        contents = "".join(f.readlines())

    with doxyfile.open("w") as f:
        f.write(replacer.replace(contents))

    if read_the_docs_build:
        # Add hook for building doxygen xml when needed
        app.connect("builder-inited", lambda _: run_doxygen(project_docs_dir))
    else:
        run_doxygen(project_docs_dir)
