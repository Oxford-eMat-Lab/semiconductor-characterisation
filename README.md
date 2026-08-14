# Semiconductor Characterisation Notebooks

Interactive Jupyter notebooks that teach the physical principles behind
semiconductor / solar-cell characterisation techniques — starting from
the basics and working through the equations, with short code examples
and plots of model data.

## Techniques

| Technique | Notebook | Open in Colab |
|---|---|---|
| Transfer Length Method (TLM) — contact resistance | [`TLM/tlm_analysis.ipynb`](TLM/tlm_analysis.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/TLM/tlm_analysis.ipynb) |
| External Quantum Efficiency (EQE) — solar cell spectral response | [`EQE/eqe_analysis.ipynb`](EQE/eqe_analysis.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/EQE/eqe_analysis.ipynb) |

## Documentation site

Key equations and background for every technique are also published as a
static site (built with MkDocs):
**https://oxford-emat-lab.github.io/semiconductor-characterisation/**

The site is built from [`docs/`](docs/) and deployed automatically by
[`.github/workflows/deploy-docs.yml`](.github/workflows/deploy-docs.yml)
on every push to `main`. To enable it: **Settings → Pages → Source:
GitHub Actions** (one-time repo setting). The workflow builds the site and
publishes it with the official Pages actions — no `gh-pages` branch is used.

## Repository structure

```
semicon_characterisation/
├── README.md
├── mkdocs.yml                  # docs site config
├── requirements.txt            # notebook dependencies
├── requirements-docs.txt       # docs site dependencies
├── docs/                       # MkDocs pages (index, tlm.md, eqe.md, assets/)
├── .github/workflows/
│   └── deploy-docs.yml         # builds + deploys the MkDocs site to GitHub Pages
├── tools/
│   ├── nb2md.py                # notebook -> MkDocs markdown converter
│   └── build_docs.sh           # regenerate the notebook-derived pages + build
├── TLM/
│   ├── tlm_analysis.ipynb      # the notebook (executed, outputs embedded)
│   ├── tlm_helper.py           # all TLM physics + fitting functions
│   ├── build_notebook.py       # source of the notebook - edit this, not the .ipynb
│   ├── make_figures.js         # builds figures.pptx from native PowerPoint shapes
│   ├── export_figures.sh       # figures.pptx -> figures/*.jpg -> docs/assets/
│   ├── figures.pptx, figures/  # the diagrams used in the notebook
│   └── contact_res_v2.ipynb    # superseded first version, kept for reference
└── EQE/
    ├── eqe_analysis.ipynb
    ├── eqe_helper.py           # all EQE physics functions used by the notebook
    ├── build_notebook.py       # source of the notebook - edit this, not the .ipynb
    ├── make_figures.js, export_figures.sh
    ├── figures.pptx, figures/
    └── literature/             # background reading used to build the notebook
```

Both techniques follow the same layout: a **helper module** holding every
physics and fitting function, a **`build_notebook.py`** that defines the
notebook cell by cell and executes it, and a **PowerPoint figure deck**
built from native shapes. The `.ipynb` files are build products — edit
`build_notebook.py` and re-run it rather than editing a notebook by hand,
or your changes will be overwritten.

## Running the notebooks locally

Two setups are described below: **macOS with `uv`** and **Windows 11 with
Miniconda**. Either one works — pick whichever matches your machine.

In both cases, launch Jupyter from the repository root. Jupyter sets each
notebook's working directory to the folder containing it, so
`import eqe_helper` / `import tlm_helper` and the `./figures/...` image
paths resolve automatically when you open either notebook.

### macOS — `uv`

[`uv`](https://docs.astral.sh/uv/) is a fast Python package/environment
manager. Install it once:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
# or, if you use Homebrew:
brew install uv
```

Then, from the repository root:

```bash
cd semicon_characterisation

uv venv                                 # create .venv (uv fetches Python if needed)
source .venv/bin/activate               # activate it
uv pip install -r requirements.txt      # numpy, scipy, matplotlib, jupyterlab, ...

jupyter lab                             # opens in your browser
```

Next time, you only need `source .venv/bin/activate` before `jupyter lab`.

To pin a specific Python version instead of the system one:

```bash
uv venv --python 3.12
```

<details>
<summary>Alternative: run without activating the environment</summary>

`uv run` resolves dependencies and runs the command in one step, so you can
skip the activate step entirely:

```bash
uv run --with-requirements requirements.txt jupyter lab
```

</details>

### Windows 11 — Miniconda

Install Miniconda using the
[Windows graphical installer](https://www.anaconda.com/docs/getting-started/miniconda/install/windows-gui-install),
then open **Anaconda Prompt** from the Start menu and run:

```bat
cd path\to\semicon_characterisation

:: create and activate the environment
conda create -n semicon python=3.12
conda activate semicon

:: install the dependencies
pip install -r requirements.txt

:: open Jupyter in your browser
jupyter lab
```

Next time, you only need `conda activate semicon` before `jupyter lab`.

If you would rather install the scientific stack through conda-forge
(instead of pip), replace the `pip install` line with:

```bat
conda install -c conda-forge numpy pandas matplotlib scipy jupyterlab
```

> **Note:** use the **Anaconda Prompt**, not the plain Windows Command
> Prompt or PowerShell — `conda activate` only works in a shell where conda
> has been initialised. (To use PowerShell instead, run `conda init powershell`
> once and restart the terminal.)

### Google Colab

No installation needed — use the Colab badges in the table above. Colab
already provides numpy, scipy, matplotlib and pandas.

Colab opens the notebook file on its own, *without* the rest of the
repository, so the first cell of each notebook clones this repo and
switches into its folder to make the helper module and `figures/`
available. Run the cells in order and it is handled automatically — but
this requires the repository to be public on GitHub.

## Building the documentation site locally

The technique pages on the site **are** the notebooks: `docs/eqe.md` and
`docs/tlm.md` are generated from the `.ipynb` files by
[`tools/nb2md.py`](tools/nb2md.py), so the published page carries the same
explanations, numbered equations, code cells and figures. Only
`docs/index.md` is hand-written.

```bash
pip install -r requirements-docs.txt
./tools/build_docs.sh serve    # regenerate pages, then serve on :8000
./tools/build_docs.sh          # regenerate pages, then mkdocs build --strict
```

Do not edit `docs/eqe.md` or `docs/tlm.md` directly — they are overwritten.
Edit the notebook instead — that means editing
`EQE/build_notebook.py` or `TLM/build_notebook.py` and re-running it — then
re-run `tools/build_docs.sh`. Notebook cells tagged `hide-in-docs` are
skipped on the site; both notebooks use this to hide their Colab setup
cell.

The GitHub Actions workflow regenerates the pages before deploying, so the
site cannot drift from the notebooks.

### Why not the mkdocs-jupyter plugin?

[`mkdocs-jupyter`](https://github.com/danielfrg/mkdocs-jupyter) renders
`.ipynb` files directly and is a perfectly good option. To use it instead,
add it to `requirements-docs.txt` and set:

```yaml
plugins:
  - mkdocs-jupyter:
      include_source: true
      execute: false          # notebooks are already executed

nav:
  - External Quantum Efficiency (EQE): EQE/eqe_analysis.ipynb
```

The converter is used here for two reasons. It produces plain markdown, so
the site does not depend on the plugin API — relevant because MkDocs 2.0 is
slated to drop plugin support entirely (which would also affect the
Material theme). And it needs nothing beyond the Python standard library,
so the docs build has no extra dependency.

## Figures

Each technique's diagrams live in a `figures.pptx`, one per slide, built
from **native PowerPoint shapes** (rectangles, arrows, text boxes) rather
than flat images — so they can be edited directly in PowerPoint.

[`EQE/figures.pptx`](EQE/figures.pptx):

| Slide | Exported as | Shows |
|---|---|---|
| 1 | `fig_cell_structure.jpg` | c-Si cell cross-section, depth axis |
| 2 | `fig_optical_losses.jpg` | photon accounting; EQE vs IQE |
| 3 | `fig_collection_efficiency.jpg` | penetration depth vs collection |
| 4 | `fig_albsf_vs_perc.jpg` | Al-BSF vs PERC rear side |
| 5 | `fig_measurement_setup.jpg` | DSR measurement schematic |

[`TLM/figures.pptx`](TLM/figures.pptx):

| Slide | Exported as | Shows |
|---|---|---|
| 1 | `fig_tlm_resistance_chain.jpg` | the series chain behind a two-probe reading |
| 2 | `fig_tlm_structure.jpg` | the TLM test pattern, top view |
| 3 | `fig_tlm_crowding.jpg` | current crowding and the transfer length |
| 4 | `fig_tlm_plot.jpg` | what each feature of the fitted line gives |
| 5 | `fig_tlm_measurement.jpg` | four-wire sweep on one contact pair |

After editing a deck, re-export the JPEGs used by the notebook and the
docs site:

```bash
cd TLM        # or cd EQE
./export_figures.sh      # needs LibreOffice + poppler-utils + ImageMagick
```

The script writes `<TECHNIQUE>/figures/*.jpg` and copies them to
`docs/assets/` (figure names are prefixed per technique so the two decks
cannot overwrite each other there). To regenerate the deck itself from
source instead, run `node make_figures.js` in that folder, then re-run the
export.

## License

Add a license (e.g. MIT) before publishing, if this repository is meant
to be open-source.
