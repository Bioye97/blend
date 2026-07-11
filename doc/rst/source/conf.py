import importlib.util

project = "BLEND"
author = "Rasheed Ajala"
copyright = "2026, Rasheed Ajala"

release = "2.0.0"
version = release

if importlib.util.find_spec("sphinx_rtd_theme") is None:
    raise RuntimeError(
        "BLEND documentation uses the Read the Docs theme. "
        "Install it with 'python -m pip install sphinx-rtd-theme' "
        "or 'conda install -c conda-forge sphinx_rtd_theme'."
    )

extensions = ["sphinx_rtd_theme"]
templates_path = ["_templates"]
exclude_patterns = []

html_theme = "sphinx_rtd_theme"
html_title = "BLEND Documentation"
html_static_path = ["_static"]
html_copy_source = False
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "style_external_links": True,
}

latex_documents = [
    ("index", "blend.tex", "BLEND Documentation", author, "manual"),
]

latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "10pt",
}
