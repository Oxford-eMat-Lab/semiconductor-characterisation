# Semiconductor Characterisation Techniques

Interactive Jupyter notebooks that explain the physical principles behind
common semiconductor / solar-cell characterisation techniques, worked
through from first principles with equations, short code examples, and
plots of model data.

Each technique has:

- a **Jupyter notebook** (runnable locally or in Google Colab) with the
  full derivation, code, and plots,
- a **helper Python module** with the underlying physics functions, kept
  separate so the notebook itself stays short and readable,
- a **page on this site**, which *is* that notebook — the same
  explanations, numbered equations, code cells and figures, rendered for
  reading without running anything.

## Techniques

| Technique | What it extracts | Notebook | Docs |
|---|---|---|---|
| Transfer Length Method (TLM) | Sheet resistance, contact resistance, transfer length, specific contact resistivity | [`TLM/tlm_analysis.ipynb`](https://github.com/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/TLM/tlm_analysis.ipynb) | [TLM](tlm.md) |
| External Quantum Efficiency (EQE) | Spectral response, absorption/collection behaviour, short-circuit current density | [`EQE/eqe_analysis.ipynb`](https://github.com/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/EQE/eqe_analysis.ipynb) | [EQE](eqe.md) |

See the repository [README](https://github.com/Oxford-eMat-Lab/semiconductor-characterisation#readme)
for Google Colab links and setup instructions.
