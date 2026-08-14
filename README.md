# Semiconductor Characterisation Notebooks

Interactive Jupyter notebooks that teach the physical principles behind
semiconductor / solar-cell characterisation techniques — starting from
the basics and working through the equations, with short code examples
and plots of model data.

> Replace `YOUR_GH_USERNAME` below with the actual GitHub username/org
> once this repository is pushed to GitHub, so the Colab and docs links
> resolve correctly.

## Techniques

| Technique | Notebook | Open in Colab |
|---|---|---|
| Transfer Length Method (TLM) — contact resistance | [`TLM/contact_res_v2.ipynb`](TLM/contact_res_v2.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_GH_USERNAME/semicon_characterisation/blob/main/TLM/contact_res_v2.ipynb) |
| External Quantum Efficiency (EQE) — solar cell spectral response | [`EQE/eqe_analysis.ipynb`](EQE/eqe_analysis.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_GH_USERNAME/semicon_characterisation/blob/main/EQE/eqe_analysis.ipynb) |

## Documentation site

Key equations and background for every technique are also published as a
static site (built with MkDocs):
**https://YOUR_GH_USERNAME.github.io/semicon_characterisation/**

The site is built from [`docs/`](docs/) and deployed automatically by
[`.github/workflows/deploy-docs.yml`](.github/workflows/deploy-docs.yml)
on every push to `main` that touches `docs/` or `mkdocs.yml`. To enable
it: **Settings → Pages → Source: GitHub Actions** (one-time repo setting).

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
├── TLM/
│   ├── contact_res_v2.ipynb
│   ├── tlm2.jpg, tlm3.jpg       # figures used in the notebook
│   └── figs.pptx                 # source slides for the figures
└── EQE/
    ├── eqe_analysis.ipynb
    ├── eqe_helper.py            # all EQE physics functions used by the notebook
    ├── figures/                  # figure placeholders (replace with real JPEGs)
    └── literature/                # background reading used to build the notebook
```

## Running the notebooks locally

Two setups are described below: **macOS with `uv`** and **Windows 11 with
Miniconda**. Either one works — pick whichever matches your machine.

In both cases, launch Jupyter from the repository root. Jupyter sets each
notebook's working directory to the folder containing it, so
`import eqe_helper` and the `./figures/...` image paths resolve
automatically when you open `EQE/eqe_analysis.ipynb`.

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
repository, so the first cell of the EQE notebook clones this repo and
switches into the `EQE/` folder to make `eqe_helper.py` and `figures/`
available. Run the cells in order and it is handled automatically — but
this only works once the repository is public on GitHub and the
`YOUR_GH_USERNAME` placeholder in that cell has been replaced.

## Building the documentation site locally

```bash
pip install -r requirements-docs.txt
mkdocs serve   # http://127.0.0.1:8000
```

## Figures

Diagram placeholders (`EQE/figures/*.jpg`) are auto-generated stand-ins
with a "TO BE REPLACED" label. Replace them with real JPEGs exported from
PowerPoint, keeping the same filenames so the notebook and docs site pick
them up automatically:

- `fig_cell_structure.jpg` — Al-BSF vs. PERC solar cell cross-section
- `fig_collection_efficiency.jpg` — collection-efficiency profile sketch
- `fig_albsf_vs_perc.jpg` — Al-BSF vs. PERC rear-side comparison
- `fig_measurement_setup.jpg` — DSR measurement setup schematic

## License

Add a license (e.g. MIT) before publishing, if this repository is meant
to be open-source.
