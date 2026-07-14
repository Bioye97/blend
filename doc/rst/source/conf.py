import importlib.util
import os
from pathlib import Path

from docutils import nodes
from docutils.parsers.rst import Directive

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
html_css_files = ["blend.css"]
html_js_files = ["blend.js"]
html_copy_source = False
html_show_sourcelink = False
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "style_external_links": True,
}
html_context = {
    "display_github": True,
    "github_user": "Bioye97",
    "github_repo": "blend",
    "github_version": "devel-2.0",
    "conf_py_path": "/doc/rst/source/",
    "blend_current_version": "devel-2.0",
    "blend_versions": [
        ("devel-2.0", "https://bioye97.github.io/blend/"),
    ],
}

latex_documents = [
    ("index", "blend.tex", "BLEND Documentation", author, "manual"),
]

latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "10pt",
    "fontpkg": "",
}


class BlendUsageDirective(Directive):
    required_arguments = 1
    optional_arguments = 0
    has_content = False

    def run(self):
        module = self.arguments[0]
        usage_dir = os.environ.get("BLEND_DOC_USAGE_DIR")
        if not usage_dir:
            warning = self.state_machine.reporter.warning(
                "BLEND_DOC_USAGE_DIR is not set; generated usage is unavailable",
                line=self.lineno,
            )
            return [warning]

        usage_path = Path(usage_dir) / f"{module}.txt"
        try:
            text = usage_path.read_text(encoding="utf-8")
        except OSError as exc:
            warning = self.state_machine.reporter.warning(
                f"Could not read generated usage file {usage_path}: {exc}",
                line=self.lineno,
            )
            return [warning]

        literal = nodes.literal_block(text, text)
        literal["language"] = "text"
        return [literal]


def setup(app):
    app.add_directive("blend-usage", BlendUsageDirective)
    return {"version": "1.0", "parallel_read_safe": True}
