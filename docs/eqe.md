# External Quantum Efficiency (EQE)

Full runnable version: [`EQE/eqe_analysis.ipynb`](https://github.com/YOUR_GH_USERNAME/semicon_characterisation/blob/main/EQE/eqe_analysis.ipynb)
([open in Colab](https://colab.research.google.com/github/YOUR_GH_USERNAME/semicon_characterisation/blob/main/EQE/eqe_analysis.ipynb)).

EQE quantifies how efficiently a solar cell converts incident photons of a
given wavelength into collected charge carriers. It is used to extract
absorption behaviour, carrier collection properties (diffusion length,
surface recombination), the short-circuit current density $J_{sc}$, and
optical losses.

<div align="center" markdown>
  <img src="../assets/fig_cell_structure.jpg" width="500">
</div>

## Photon energy and absorption

$$
E_{\text{phot}}(\lambda) = \frac{hc}{\lambda}
$$
Silicon's indirect bandgap ($E_g = 1.12\ \text{eV}$) corresponds to
$\lambda \approx 1107\ \text{nm}$, the absorption edge. The relevant range
for c-Si solar cells is $\approx 300$-$1200\ \text{nm}$.

Light intensity follows the Lambert-Beer law:
$$
\Phi(\lambda, z) = \Phi_0(\lambda)\, e^{-\alpha(\lambda) z}, \qquad
L_\alpha(\lambda) = \frac{1}{\alpha(\lambda)}
$$
$\alpha(\lambda)$ (the absorption coefficient) spans six orders of
magnitude across the useful range, because silicon's indirect gap makes
absorption near $E_g$ phonon-assisted and gradual.

## Generation and collection

$$
g(\lambda, z) = (1-R)\,\Phi_0(\lambda)\,\alpha(\lambda)\, e^{-\alpha(\lambda) z}
$$
Collection efficiency in the base region (thickness $W$, diffusion length
$L$, rear recombination velocity $S$):
$$
\eta_c(z) = \cosh\!\left(\frac{z}{L}\right) - \frac{L}{L_{\text{eff}}}\sinh\!\left(\frac{z}{L}\right),
\qquad
L_{\text{eff}} = L\,\frac{S\sinh(W/L) + D\cosh(W/L)}{S\cosh(W/L) + D\sinh(W/L)}
$$

<div align="center" markdown>
  <img src="../assets/fig_collection_efficiency.jpg" width="450">
</div>

## EQE and short-circuit current

$$
\text{EQE}(\lambda) = \int_0^{W} g(\lambda,z)\,\eta_c(z)\, dz \Big/ \Phi_0(\lambda)
= T_{\text{ext}}(\lambda)\cdot \text{IQE}(\lambda)
$$
$$
j_{sc} = q \int \Phi_0(\lambda)\, \text{EQE}(\lambda)\, d\lambda
$$

**Short wavelengths** (300-500 nm) probe the front surface / ARC / emitter
(blue response). **Long wavelengths** (900-1200 nm) probe the bulk
diffusion length and rear surface (red / near-IR response) — this is why a
PERC cell (better rear passivation) typically shows higher near-IR EQE
than an Al-BSF cell.

<div align="center" markdown>
  <img src="../assets/fig_albsf_vs_perc.jpg" width="500">
</div>

## Measuring EQE: spectral responsivity

EQE is not measured directly. The spectral responsivity
$s(\lambda) = \text{EQE}(\lambda)\cdot q\lambda / hc$ (A/W) is measured
instead. A **differential** measurement is used: steady white bias light
sets realistic injection conditions, and a small chopped monochromatic
signal gives the differential spectral responsivity (DSR),
$$
\tilde s(\lambda, E_{\text{bias}}) = \frac{\Delta j_{sc}(\lambda)}{\Delta E_\lambda(\lambda)}
$$
calibrated against a reference cell of known DSR.

<div align="center" markdown>
  <img src="../assets/fig_measurement_setup.jpg" width="500">
</div>

## Absolute scaling

The EQE from a DSR measurement is only relative (an unknown calibration
ratio remains). It is scaled to match an independently measured $J_{sc}$
(e.g. sun simulator):
$$
f_{sc} = \frac{j_{sc,\text{exp}}}{j_{sc,\text{calc}}}, \qquad
\text{EQE}_{\text{abs}} = f_{sc}\cdot \text{EQE}_{\text{rel}}
$$

## Linearity / bias-ramp measurements

For a linear cell, $\tilde s$ is independent of bias irradiance
$E_{\text{bias}}$. IEC 60904-8 defines simplified single-bias-point
procedures (e.g. $\approx 300\ \text{W/m}^2$) that approximate the full
integral
$s_{\text{STC}}(\lambda) = \int_0^{1000\,\text{W/m}^2} \tilde s\, dE_{\text{bias}}$
well for typical silicon cells; deviations are mostly below 10%, largest
around 1000 nm for nonlinear front-junction cells (Bothe et al., 2018).

## Assumptions and limitations

- one-dimensional carrier transport, uniform doping,
- simplified analytic reflectance and absorption-coefficient models,
- an analytic approximation of the AM1.5G spectrum (not tabulated
  IEC 60904-3 data),
- no series-resistance effects beyond the toy bias-ramp model.

## References

1. C. Schinke, S. Schädlich, T. Gewohn, D. Hinken, *Analysis of the Quantum
   Efficiency of Silicon Solar Cells*, Leibniz Universität Hannover / ISFH, 2019.
2. K. Bothe, D. Hinken, B. Min, C. Schinke, "Accuracy of Simplifications
   for Spectral Responsivity Measurements of Solar Cells," *IEEE J.
   Photovolt.* 8(2), 611-620, 2018.
3. Quokka3 Modelling Guide — optical modelling (EQE = T_ext × IQE).
4. IEC 60904-8, *Measurement of spectral responsivity of a PV device*.
