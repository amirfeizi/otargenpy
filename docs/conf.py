"""Sphinx configuration for otargenpy documentation."""

project = "otargenpy"
copyright = "2024-2026, Amir Feizi"
author = "Amir Feizi"
release = "2.0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstrings = True
napoleon_numpy_docstrings = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = "otargenpy"
