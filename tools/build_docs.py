#!/usr/bin/env python3
"""
Regenerates the notebook-derived pages of the docs site, then builds it.

    python tools/build_docs.py            # regenerate pages + mkdocs build
    python tools/build_docs.py serve      # regenerate pages + mkdocs serve

Pure Python and cross-platform, so it needs no bash. `tools/build_docs.sh`
is a thin wrapper around this file, which keeps CI (and anyone on a shell)
working without a second copy of the logic to keep in step.

docs/index.md is hand-written; docs/eqe.md, docs/tlm.md and docs/kpspv.md
are generated from the notebooks and must not be edited directly - edit the
notebook and run this again.

Adding a technique is one entry in TECHNIQUES below.
"""
import os
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

from nb2md import convert  # noqa: E402  (needs the path set above)

REPO_URL = os.environ.get(
    "REPO_URL",
    "https://github.com/Oxford-eMat-Lab/semiconductor-characterisation",
)

# folder, notebook stem, docs page stem
TECHNIQUES = [
    ("EQE",   "eqe_analysis",   "eqe"),
    ("TLM",   "tlm_analysis",   "tlm"),
    ("KPSPV", "kpspv_analysis", "kpspv"),
]

ASSETS_DIR = "docs/assets/nb"
ASSETS_URL = "assets/nb"


def regenerate():
    # Paths passed to convert() must be RELATIVE to the repo root: nb2md
    # writes the notebook path into the page header and into the Colab and
    # GitHub links, so an absolute path leaks the build machine's directory
    # layout into the published page and breaks both links.
    for folder, nb_stem, page in TECHNIQUES:
        notebook = f"{folder}/{nb_stem}.ipynb"
        if not (ROOT / notebook).exists():
            sys.exit(f"missing notebook: {notebook}")
        helper = f"{page}_helper.py"
        convert(
            notebook,
            f"docs/{page}.md",
            ASSETS_DIR,
            ASSETS_URL,
            [
                ("./figures/", "../assets/"),
                (f"./{helper}", f"{REPO_URL}/blob/main/{folder}/{helper}"),
            ],
            REPO_URL,
        )


def run_mkdocs(command):
    """Invoke mkdocs in-process, so it uses this interpreter's environment."""
    argv = ["mkdocs", command]
    if command == "build":
        argv.append("--strict")
    sys.argv = argv
    try:
        runpy.run_module("mkdocs", run_name="__main__")
    except SystemExit as exc:
        raise SystemExit(exc.code)
    except ImportError:
        sys.exit(
            "mkdocs is not installed in this Python environment.\n"
            f"  {sys.executable} -m pip install -r requirements-docs.txt"
        )


def main(argv):
    command = argv[1] if len(argv) > 1 else "build"
    if command not in ("build", "serve"):
        sys.exit(f"usage: {Path(__file__).name} [build|serve]")
    os.chdir(ROOT)
    regenerate()
    run_mkdocs(command)


if __name__ == "__main__":
    main(sys.argv)
