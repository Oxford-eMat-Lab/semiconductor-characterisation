<!-- GENERATED FILE - do not edit. Produced from CGV/cgv_analysis.ipynb by tools/nb2md.py (see tools/build_docs.sh). -->

!!! info "Generated from a Jupyter notebook"
    This page is `CGV/cgv_analysis.ipynb`, rendered with its stored outputs.
    [Run it in Google Colab](https://colab.research.google.com/github/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/CGV/cgv_analysis.ipynb) or
    [view the notebook on GitHub](https://github.com/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/CGV/cgv_analysis.ipynb).


# Capacitance-Voltage and Conductance-Voltage (C-V/G-V)

A C-V/G-V measurement puts a small AC signal on top of a DC bias and reads
the resulting admittance. The capacitance traces how far the depletion
layer under a gate moves as bias changes, so its slope gives a doping
profile. The conductance traces energy lost to charge hopping in and out
of interface states, so its frequency dependence gives an interface trap
density. One sweep, two very different stories, because one component is
reactive and the other is dissipative.

This notebook builds a MOS (metal-oxide-semiconductor) capacitor from
first principles, uses its C-V curve to extract a doping profile, shows
that the same algebra collapses to the classic Mott-Schottky line the
moment the oxide is removed, and then uses the G-V data from the same
sweep to find an interface trap density.

| Section | Question answered |
|---|---|
| 1 | What do we actually do when we sweep bias and read capacitance? |
| 2 | What is a MOS capacitor, and why does it behave like two capacitors in series? |
| 3 | What happens at the semiconductor surface as the gate bias changes? |
| 4 | Why is the semiconductor's share of the capacitance itself bias-dependent? |
| 5 | Why does measurement frequency change the shape of the curve? |
| 6 | What does the ideal C-V curve look like, feature by feature? |
| 7 | Where does the voltage axis sit - flat-band voltage and oxide charge |
| 8 | Where does inversion begin - threshold voltage |
| 9 | How do we turn C(V) into a doping profile? |
| 10 | What happens without an oxide? Schottky C-V and the Mott-Schottky line |
| 11 | What is the conductance method actually measuring? |
| 12 | How do we get the semiconductor's true parallel conductance out of what the meter reports? |
| 13 | What does the conductance peak tell us about interface traps? |
| 14 | What ruins a real measurement, and how do you know? |
| 15 | Assumptions and limitations |

<div align="center">
   <img src="../assets/fig_cgv_mos_structure.jpg" width="720">
</div>

Equations are numbered (1), (2), ... and referred to by those numbers
throughout. All physics functions live in
[`https://github.com/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/CGV/cgv_helper.py`](https://github.com/Oxford-eMat-Lab/semiconductor-characterisation/blob/main/CGV/cgv_helper.py), so the notebook itself stays short;
that module's docstrings point back to these equation numbers.

```python
import numpy as np
import matplotlib.pyplot as plt
import cgv_helper as cgv

plt.rcParams.update({'font.size': 12})
```

## 1. What do we actually do when we sweep bias and read capacitance?

You apply a DC voltage to a gate electrode, add a small AC signal on top
of it (a few tens of mV), and an LCR meter reports the in-phase and
out-of-phase current at that bias: a capacitance and a conductance. Sweep
the bias and you get $C(V)$ and $G(V)$.

That is the entire experiment. Every quantity this notebook extracts -
oxide thickness, doping profile, flat-band voltage, threshold voltage,
interface trap density - is *inferred* from those two curves. Nothing
else is measured directly, and that is worth stating plainly: this is the
C-V/G-V version of the measured-vs-inferred distinction that runs through
every notebook in this repository. EQE infers efficiency from spectral
response; TLM infers contact resistance from an intercept; here, an
admittance is measured and a doping profile, a trap density, and half a
dozen voltages in between are all inferred from it through a model.

## 2. What is a MOS capacitor, and why does it behave like two capacitors in series?

A MOS capacitor is a metal gate, a thin oxide, and a doped silicon
substrate, with an ohmic contact on the back. The gate and the substrate
are the two plates; the oxide is the dielectric between them - except
that one of the two plates is a semiconductor, and a semiconductor's
surface charge moves in response to the field the way a metal's never
does.

The oxide capacitance is fixed by geometry alone:

$$
C_{\text{ox}} = \frac{\varepsilon_{\text{ox}} A}{t_{\text{ox}}} \tag{1}
$$

This never changes with bias. Every feature in a C-V curve comes from the
*other* plate: the semiconductor's depletion layer growing and shrinking
underneath the gate, which changes the effective thickness of the
dielectric between the gate and the bulk silicon. The two capacitances -
fixed oxide, bias-dependent semiconductor - sit in series, and Eq. (1) is
the one number on the whole curve you can check by inspection: it should
match the flattest part of the plot in strong accumulation (Sec. 6).

```python
tox_nm = 10.0        # oxide thickness
area_cm2 = 1e-4       # 100 um x 100 um gate pad
Cox_pa = cgv.oxide_capacitance_per_area(tox_nm)
Cox = cgv.oxide_capacitance(area_cm2, tox_nm)
print(f"C_ox (per area) = {Cox_pa*1e6:.2f} uF/cm^2")
print(f"C_ox (this device, A={area_cm2*1e8:.0f} um^2) = {Cox*1e12:.2f} pF")
```

```text
C_ox (per area) = 0.35 uF/cm^2
C_ox (this device, A=10000 um^2) = 34.53 pF
```

## 3. What happens at the semiconductor surface as the gate bias changes?

Define the surface potential $\phi_s$: the band bending at the
semiconductor-oxide interface, relative to the neutral bulk, with
$\phi_s=0$ at flat band. This notebook works a p-type substrate as its
running example (n-type is a mirror image, and every function below takes
a `dopant_type` argument).

**Sign convention - declared once, used everywhere below.** For a p-type
substrate:

- $\phi_s < 0$ is **accumulation**: majority holes pulled to the surface.
- $\phi_s > 0$ is **depletion**: holes pushed away, exposing ionised
  acceptors.
- $\phi_s > 2\phi_F$ is **(strong) inversion**: enough band bending that
  the surface behaves as if it were n-type, where the bulk Fermi
  potential is

$$
\phi_F = V_t \ln\!\left(\frac{N_A}{n_i}\right) \tag{2}
$$

positive for p-type. Increasing the gate bias moves a p-type surface from
accumulation toward inversion; for n-type it is the reverse, and $\phi_F$
itself is negative.

This is the same band-bending picture the Kelvin probe / surface
photovoltage notebook uses - there, $\phi_s$ is inferred from a null
voltage or from illumination; here, it is driven directly by an applied
gate bias. The electrostatics underneath is identical.

## 4. Why is the semiconductor's share of the capacitance itself bias-dependent?

Solving Poisson's equation with Boltzmann statistics gives the exact
semiconductor surface field, valid all the way from accumulation through
depletion into inversion (not just the depletion approximation) [[5]](#ref5):

$$
E_s = 2\,\text{sign}(\phi_s)\sqrt{\frac{qn_iV_t}{\varepsilon_s}}
\sqrt{\cosh\!\left(\frac{\phi_s-\phi_F}{V_t}\right)
+ \frac{\phi_s}{V_t}\sinh\!\left(\frac{\phi_F}{V_t}\right)
- \cosh\!\left(\frac{\phi_F}{V_t}\right)}
$$

and, from Gauss's law, the space-charge density per unit area:

$$
Q_{sc}(\phi_s) = -\varepsilon_s E_s \tag{3}
$$

$Q_{sc}$ is a *decreasing* function of $\phi_s$: more band bending always
drives the semiconductor's own charge more negative (positive in
accumulation, negative in depletion, for p-type). The capacitance that
matters is the gate's response, $-Q_{sc}$, so

$$
C_s(\phi_s) = -\frac{dQ_{sc}}{d\phi_s} \tag{4}
$$

a **differential** quantity - the slope of Eq. (3), not the charge
divided by the potential. Two useful checks: at flat band, $C_s$ reduces
to $\varepsilon_s/L_D$, with the Debye length

$$
L_D = \sqrt{\frac{\varepsilon_s V_t}{qN}} \tag{5}
$$

and, well inside depletion (small $\phi_s$, far from inversion), $Q_{sc}$
matches the familiar depletion-approximation width

$$
W(\phi_s) \approx \sqrt{\frac{2\varepsilon_s\phi_s}{qN}} \tag{6}
$$

to a few percent. Both checks are in `cgv_helper_checks.py`.

```python
Na = 1e16   # p-type substrate, cm^-3
phi_s = np.linspace(-0.4, 1.0, 400)
Qsc = cgv.space_charge_density(phi_s, Na, 'p')
Cs = cgv.semiconductor_capacitance(phi_s, Na, 'p')
phiF = cgv.fermi_potential(Na, 'p')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.plot(phi_s, Qsc*1e6)
ax1.axvline(0, color='gray', lw=0.8)
ax1.axvline(2*phiF, color='gray', ls='--', lw=0.8)
ax1.set_xlabel(r'surface potential $\phi_s$ (V)')
ax1.set_ylabel(r'$Q_{sc}$ ($\mu$C/cm$^2$)')
ax1.set_title('space charge (Eq. 3)')

ax2.semilogy(phi_s, Cs*1e6)
ax2.axvline(0, color='gray', lw=0.8, label='flat band')
ax2.axvline(2*phiF, color='gray', ls='--', lw=0.8, label=r'$2\phi_F$ (inversion onset)')
ax2.set_xlabel(r'surface potential $\phi_s$ (V)')
ax2.set_ylabel(r'$C_s$ ($\mu$F/cm$^2$)')
ax2.set_title('semiconductor capacitance (Eq. 4)')
ax2.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 1](assets/nb/cgv_analysis_01.png)

$C_s$ spans nearly five decades between accumulation and depletion.
Accumulation piles up majority carriers right at the interface, so a tiny
change in $\phi_s$ moves a lot of charge - $C_s$ is huge, and the series
combination in Sec. 6 will be dominated by $C_{\text{ox}}$ alone. In
depletion, $C_s$ falls to the value set by the depletion width, which is
why the total capacitance drops as the gate is biased away from
accumulation. Past $2\phi_F$, $C_s$ turns back up as the exponentially
growing minority-carrier (inversion) charge takes over from the slowly
saturating depletion charge - and whether that turn-up shows up in a real
measurement is entirely a question of frequency, which is Sec. 5.

## 5. Why does measurement frequency change the shape of the curve?

Majority carriers respond to the AC signal in about a dielectric
relaxation time - picoseconds, effectively instantaneous. Minority
carriers in the depletion region do not exist in equilibrium; they have
to be *generated*, thermally, one electron-hole pair at a time, and that
takes a generation lifetime - microseconds to milliseconds in a
teaching-grade device. That mismatch is the whole reason C-V curves are
measured at more than one frequency:

- **High frequency (HF).** The AC signal is far too fast for minority
  carriers to generate in step with it. Past the onset of inversion, the
  inversion charge is frozen at whatever the DC bias has managed to
  accumulate over the whole sweep, and $C_s$ is clamped at its
  depletion-limited value.
- **Low frequency / quasi-static (LF).** The AC signal (or the DC ramp
  rate, if a true DC method is used instead) is slow enough that minority
  carriers keep up. The equilibrium $C_s(\phi_s)$ of Eq. (4) - inversion
  turn-up included - is measured directly.
- **Deep depletion.** The *bias sweep itself* is faster than the
  generation lifetime, so the surface never reaches equilibrium at all:
  depletion keeps growing past where it would have stopped, and $C$ falls
  below even the HF minimum.

Section 6 builds all three from the same underlying physics, differing
only in what is allowed to happen to the minority-carrier charge.

## 6. What does the ideal C-V curve look like, feature by feature?

<div align="center">
   <img src="../assets/fig_cgv_curve_shapes.jpg" width="680">
</div>

The oxide and the semiconductor combine in series:

$$
C(V) = \frac{C_{\text{ox}}\,C_s(\phi_s(V))}{C_{\text{ox}} + C_s(\phi_s(V))}
\tag{7}
$$

with $\phi_s(V)$ found from charge balance across the structure - the
gate voltage splits between band bending and the voltage dropped across
the oxide:

$$
V_G = V_{FB} + \phi_s - \frac{Q_{sc}(\phi_s)}{C_{\text{ox}}} \tag{8}
$$

solved numerically for $\phi_s$ at each $V_G$ (`surface_potential_from_bias`).
$V_{FB}$, the flat-band voltage, is the subject of Sec. 7; take it as a
given offset for now.

```python
VFB = -0.9   # flat-band voltage, worked example
Vg = np.linspace(VFB - 3.0, VFB + 3.0, 300)

C_hf = cgv.hf_cv_curve(Vg, area_cm2, tox_nm, Na, 'p', VFB)
C_lf = cgv.lf_cv_curve(Vg, area_cm2, tox_nm, Na, 'p', VFB)
C_dd = cgv.deep_depletion_cv(Vg, area_cm2, tox_nm, Na, 'p', VFB)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(Vg, C_hf*1e12, label='high frequency', lw=2)
ax.plot(Vg, C_lf*1e12, label='low frequency / quasi-static', lw=2)
ax.plot(Vg, C_dd*1e12, label='deep depletion (swept too fast)', ls='--')
ax.axhline(Cox*1e12, color='gray', lw=0.8, label=r'$C_{ox}$')
ax.axvline(VFB, color='gray', ls=':', lw=0.8, label=r'$V_{FB}$')
ax.set_xlabel(r'gate voltage $V_G$ (V)')
ax.set_ylabel('capacitance (pF)')
ax.set_title('Eqs. (7)-(8): p-type MOS capacitor, $N_A=10^{16}$ cm$^{-3}$, 10 nm oxide')
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 2](assets/nb/cgv_analysis_02.png)

Three plateaus and a knee, read left to right. In strong
accumulation, all three curves sit at $C_{\text{ox}}$ (slightly below it,
because a 10 nm oxide is thin enough that accumulation does not fully
saturate - a real, literature-documented effect, not noise). Through
depletion, all three fall together, because the physics is identical
until minority carriers matter. Past the knee, they separate: HF freezes
at a minimum set by the maximum depletion width; LF rises back toward
$C_{\text{ox}}$ as the inversion layer forms in step with the sweep; and
deep depletion overshoots past the HF minimum because the surface never
gets the chance to invert at all - it is not swept too fast in *time*,
it is swept too fast relative to how slowly minority carriers generate.
Comparing the HF and LF curves at fixed bias is itself a diagnostic: any
difference beyond depletion is minority-carrier response, and its size
is a rough clock on the generation lifetime.

## 7. Where does the voltage axis sit - flat-band voltage and oxide charge?

$$
V_{FB} = \phi_{ms} - \frac{Q_{\text{eff}}}{C_{\text{ox}}} \tag{9}
$$

the metal-semiconductor work function difference $\phi_{ms}$, minus the
effective oxide charge divided by the oxide capacitance [[2]](#ref2). If
the oxide carried no charge at all, $V_{FB}$ would just be $\phi_{ms}$;
any additional shift is oxide charge, and Eq. (9) rearranges to read that
charge back off a measured curve:

$$
Q_{\text{eff}} = C_{\text{ox}}\left(\phi_{ms} - V_{FB}\right) \tag{10}
$$

**This is the section's own measured-vs-inferred moment**: $V_{FB}$ is
read directly off a curve, $Q_{\text{eff}}$ is inferred from how far that
curve has moved relative to the charge-free ideal.

To find $V_{FB}$ on a measured curve, compute the *ideal* flat-band
capacitance from the Debye length at $\phi_s=0$ - $C_s(0)=\varepsilon_s/L_D$
combined in series with $C_{\text{ox}}$ - and read off the bias where the
measured HF curve crosses it.

```python
Cfb_pa = cgv.flatband_capacitance_per_area(Cox_pa, Na)
Cfb = Cfb_pa * area_cm2
i_cross = np.argmin(np.abs(C_hf - Cfb))
VFB_read = Vg[i_cross]
print(f"C_FB = {Cfb*1e12:.3f} pF")
print(f"V_FB read off the curve = {VFB_read:.3f} V  (true value: {VFB:.3f} V)")

phi_ms = -0.85
Qeff = cgv.effective_oxide_charge(phi_ms, VFB_read, Cox_pa)
Neff = Qeff / cgv.Q
print(f"phi_ms (assumed) = {phi_ms:.3f} eV")
print(f"Q_eff = {Qeff*1e9:.3f} nC/cm^2  ->  N_eff = {Neff:.3e} cm^-2")
```

```text
C_FB = 14.614 pF
V_FB read off the curve = -0.890 V  (true value: -0.900 V)
phi_ms (assumed) = -0.850 eV
Q_eff = 13.801 nC/cm^2  ->  N_eff = 8.614e+10 cm^-2
```

## 8. Where does inversion begin - threshold voltage?

Threshold is defined as the bias at which $\phi_s$ reaches $2\phi_F$ -
the conventional onset of strong inversion:

$$
V_T = V_{FB} + s\left[\frac{\sqrt{4\varepsilon_s q N|\phi_F|}}{C_{\text{ox}}}\right]
+ s\cdot 2|\phi_F| \tag{11}
$$

with $s=+1$ for p-type ($V_T>V_{FB}$, more depletion needed to invert)
and $s=-1$ for n-type ($V_T<V_{FB}$) [[6]](#ref6). On the HF curve, $V_T$
sits close to (though not exactly at) the capacitance minimum - the true
minimum is set by where the exact $C_s(\phi_s)$ of Eq. (4) actually turns
over, which is what Sec. 4's checks locate numerically rather than assume
analytically.

```python
VT = cgv.threshold_voltage(Cox_pa, Na, 'p', VFB)
print(f"V_T = {VT:.3f} V   (V_FB = {VFB:.3f} V)")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(Vg, C_hf*1e12, lw=2)
ax.axvline(VFB, color='gray', ls=':', label=r'$V_{FB}$')
ax.axvline(VT, color='crimson', ls='--', label=r'$V_T$')
ax.set_xlabel(r'gate voltage $V_G$ (V)')
ax.set_ylabel('HF capacitance (pF)')
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 3](assets/nb/cgv_analysis_03.png)

```text
V_T = -0.043 V   (V_FB = -0.900 V)
```

## 9. How do we turn C(V) into a doping profile?

<div align="center">
   <img src="../assets/fig_cgv_doping_profile.jpg" width="680">
</div>

Differentiate the series combination and the algebra collapses to a
single working formula [[1]](#ref1):

$$
N(W) = \frac{2}{q\varepsilon_s A^2\left|d(1/C^2)/dV\right|} \tag{12}
$$

with a depth scale that comes straight from the measured capacitance
alone:

$$
W = \frac{\varepsilon_s A}{C} \tag{13}
$$

The absolute value in Eq. (12) is deliberate, and worth dwelling on: the
textbook form of this formula is usually written for a bias convention
where increasing voltage means increasing *reverse* bias, so $1/C^2$
falls as $V$ rises. The convention used throughout this notebook has
increasing $V_G$ deepen depletion directly (for p-type), the opposite
sign - so $1/C^2$ *rises* with $V_G$ here. $N(W)$ itself can never be
negative, so working with $|d(1/C^2)/dV|$ sidesteps having to track which
convention is in play. Get this wrong and the formula still runs; it just
hands back a negative doping density, which is a very easy mistake to
miss if you are not looking for it.

This is the notebook's TLM-equivalent section: a slope carries the
physics, extracted point by point from noisy data, and a naive finite
difference is a lot noisier than the underlying curve suggests. Compare a
raw finite difference against a smoothed one on the same data.

```python
N_true_uniform = 3e16
VFB2 = -0.8
Vg2 = np.linspace(VFB2 + 0.1, VFB2 + 1.4, 120)
C2 = cgv.synthetic_hf_cv(Vg2, area_cm2, tox_nm, N_true_uniform, 'p', VFB2, noise_frac=0.01, seed=1)

W_raw, N_raw = cgv.doping_profile_from_cv(Vg2, C2, area_cm2, smooth=False)
W_sm, N_sm = cgv.doping_profile_from_cv(Vg2, C2, area_cm2, smooth=True)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(W_raw*1e7, N_raw, '.', ms=3, alpha=0.5, label='raw finite difference')
ax.plot(W_sm*1e7, N_sm, lw=2, label='smoothed (Savitzky-Golay)')
ax.axhline(N_true_uniform, color='gray', ls='--', label='true doping')
ax.set_xlabel('depth W (nm)')
ax.set_ylabel(r'N (cm$^{-3}$)')
ax.set_title('Eqs. (12)-(13): recovered doping, 1% capacitance noise')
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 4](assets/nb/cgv_analysis_04.png)

The smoothed curve tracks the true, uniform doping level closely;
the raw finite difference of the same data swings by a factor of a few
around it, because differentiating amplifies noise. Neither curve is
"wrong" - they differentiate the same measurement - but only one is
usable for anything quantitative. This is the same lesson TLM teaches
with its fitted intercept: the transformation between what is measured
and what is reported is not free, and a doping profile inherits every bit
of noise in the raw capacitance, amplified by a derivative.

A second, sharper test: a profile that genuinely is not uniform.

```python
def two_layer_doping(V):
    # a step in doping partway through the sweep - an implanted layer over
    # a lower-doped substrate, the kind of profile that actually needs Eq. (12)
    N1, N2 = 8e16, 2e16
    return np.where(V < VFB2 + 0.5, N1, N2)

# Build synthetic C(V) for a *non-uniform* profile directly from the
# depletion-approximation relation W(N) inverted numerically at each step -
# illustrative, not exact, since Eqs. (7)-(8) assume uniform N throughout.
V_nu = np.linspace(VFB2 + 0.08, VFB2 + 1.3, 150)
N_local = two_layer_doping(V_nu)
# integrate W(V) = sqrt(2 eps_s phi_s /(q N(phi_s))) approximately by
# treating each step as locally uniform and stepping phi_s ~ (V-VFB2)
phi_s_approx = np.clip(V_nu - VFB2, 1e-4, None)
W_nu = np.sqrt(2*cgv.EPS_SI*phi_s_approx/(cgv.Q*N_local))
Cs_nu = cgv.EPS_SI / W_nu
C_nu = cgv.total_capacitance_series_per_area(Cox_pa, Cs_nu) * area_cm2
C_nu *= (1 + np.random.default_rng(2).normal(0, 0.01, size=C_nu.shape))

W_nu_ex, N_nu_ex = cgv.doping_profile_from_cv(V_nu, C_nu, area_cm2, smooth=True)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(W_nu_ex*1e7, N_nu_ex, lw=2, label='recovered')
ax.step(W_nu*1e7, N_local, where='mid', color='gray', ls='--', label='true (step profile)')
ax.set_xlabel('depth W (nm)')
ax.set_ylabel(r'N (cm$^{-3}$)')
ax.set_title('a genuinely non-uniform profile: implant over substrate')
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 5](assets/nb/cgv_analysis_05.png)

The recovered profile shows the step, but rounded off over a few
tens of nanometres rather than sharp - a real, physical effect (carriers
diffuse across the boundary and the field redistributes over a few Debye
lengths, Eq. (5)), not just a numerical artefact of the smoothing. Barnes
[[1]](#ref1) reports exactly this rounding in real GaAs profiling data.
The spatial resolution of C-V doping profiling is fundamentally limited
to a few Debye lengths, no matter how clean the measurement.

## 10. What happens without an oxide? Schottky C-V and the Mott-Schottky line

Set $C_{\text{ox}}\to\infty$ in Eqs. (7) and (12) - no series oxide
element at all, a bare metal-semiconductor (Schottky) junction - and the
general result collapses to the classic **Mott-Schottky** line:

$$
\frac{1}{C^2} = \frac{2(V_{bi}-V)}{q\varepsilon_s A^2 N} \tag{14}
$$

linear in $V$: the slope gives $N$, the intercept on the voltage axis
gives the built-in potential $V_{bi}$ [[1]](#ref1). This is deliberately
the same payoff KPSPV gives when SPV falls out of the CPD machinery as a
special case - a general result specialising to a simpler, more famous
one. It is also, structurally, the same fit TLM makes: a straight line
whose slope is well-constrained by the data and whose intercept is a more
fragile extrapolation.

The important check is that this is an *algebraic identity*, not a
numerical coincidence - the same helper code that builds the MOS curves
above, with `Cox_pa` set enormous, reproduces Eq. (14) exactly
(`cgv_helper_checks.py`, check 10).

```python
N_schottky = 3e16
Vbi_true = 0.75
V_ms = np.linspace(-8, -0.2, 60)
C_ms = cgv.synthetic_mott_schottky(V_ms, N_schottky, Vbi_true, area_cm2, noise_frac=0.03, seed=3)

N_fit, N_fit_err, Vbi_fit, Vbi_fit_err = cgv.mott_schottky_fit(V_ms, C_ms, area_cm2)
print(f"N (true)   = {N_schottky:.3e} cm^-3")
print(f"N (fit)    = {N_fit:.3e} +/- {N_fit_err:.1e} cm^-3  ({100*N_fit_err/N_fit:.1f}%)")
print(f"Vbi (true) = {Vbi_true:.3f} V")
print(f"Vbi (fit)  = {Vbi_fit:.3f} +/- {Vbi_fit_err:.3f} V  ({100*Vbi_fit_err/Vbi_true:.1f}%)")

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(V_ms, 1/C_ms**2, 'o', ms=4, alpha=0.6, label='synthetic data (3% noise)')
fit_line = 2*(Vbi_fit - V_ms)/(cgv.Q*cgv.EPS_SI*area_cm2**2*N_fit)
ax.plot(V_ms, fit_line, lw=2, label='linear fit')
ax.set_xlabel('bias V (V)')
ax.set_ylabel(r'$1/C^2$ (F$^{-2}$)')
ax.set_title('Eq. (14): Mott-Schottky analysis')
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 6](assets/nb/cgv_analysis_06.png)

```text
N (true)   = 3.000e+16 cm^-3
N (fit)    = 3.013e+16 +/- 3.3e+14 cm^-3  (1.1%)
Vbi (true) = 0.750 V
Vbi (fit)  = 0.767 +/- 0.053 V  (7.1%)
```

$N$ comes back to within a percent or two of the noise level; $V_{bi}$,
extrapolated out to where the line crosses zero, carries a noticeably
larger relative uncertainty for the same input noise - the slope is
constrained by every point on the line, the intercept only by how far you
have to extrapolate to reach it. Exactly the TLM lesson, in a different
notebook.

## 11. What is the conductance method actually measuring?

<div align="center">
   <img src="../assets/fig_cgv_conductance_method.jpg" width="680">
</div>

Interface states sit at energy levels inside the silicon band gap, right
at the oxide interface. At a given bias, states near the Fermi level
capture and emit carriers on some characteristic timescale $\tau_{it}$.
When the AC signal's angular frequency $\omega$ is comparable to
$1/\tau_{it}$, the traps cannot quite keep up: their charging lags the
signal, and that lag dissipates energy. A genuine loss mechanism - it
shows up in $G$, the out-of-phase component, not in $C$.

The measured admittance at any one bias and frequency is

$$
Y_m = G_m + j\omega C_m \tag{15}
$$

[[3]](#ref3). Sweep frequency at fixed bias, and the interface-trap
signature is a peak in $G_m/\omega$ - but only *after* removing the
series oxide capacitance from the circuit, which is Sec. 12.

## 12. How do we get the semiconductor's true parallel conductance out of what the meter reports?

<div align="center">
   <img src="../assets/fig_cgv_equivalent_circuit.jpg" width="680">
</div>

$C_{\text{ox}}$ still sits in series with the semiconductor branch in the
equivalent circuit, so the raw $G_m/\omega$ is not yet the quantity that
peaks cleanly with interface-trap physics. Converting the measured
admittance through the known $C_{\text{ox}}$ gives the semiconductor
branch's true parallel conductance [[3]](#ref3):

$$
\frac{G_p}{\omega} = \frac{\omega C_{\text{ox}}^2 G_m}
{G_m^2 + \omega^2(C_{\text{ox}}-C_m)^2} \tag{16}
$$

Skipping this step, or using an inaccurate $C_{\text{ox}}$, is the single
biggest source of a wrong $D_{it}$ - the error in $C_{\text{ox}}$ propagates
directly into Eq. (16) and from there directly into everything Sec. 13
extracts.

```python
Dit_true = 5e11    # cm^-2 (areal density, single-level teaching model)
tau_it_true = 2e-5  # s
freqs = np.logspace(2, 6, 100)
Cm, Gm = cgv.synthetic_conductance_sweep(
    freqs, Dit_true, tau_it_true, area_cm2, tox_nm, noise_frac=0.02, seed=4)
omega = 2*np.pi*freqs
Gp_over_w = cgv.admittance_to_parallel(Cm, Gm, omega, Cox)

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogx(freqs, Gm/omega*1e12, 'o-', ms=3, alpha=0.5, label=r'raw $G_m/\omega$ (uncorrected)')
ax.semilogx(freqs, Gp_over_w*1e12, 'o-', ms=3, label=r'$G_p/\omega$ (Eq. 16, corrected)')
ax.set_xlabel('frequency (Hz)')
ax.set_ylabel(r'$G/\omega$ (pF)')
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 7](assets/nb/cgv_analysis_07.png)

The corrected curve peaks at a lower frequency and a higher value
than the raw one - the series $C_{\text{ox}}$ shortens the semiconductor
branch's effective time constant when left uncorrected, exactly as
Sec. 12 warns. Only the corrected curve is the one Sec. 13 fits.

## 13. What does the conductance peak tell us about interface traps?

The simplest possible model - a single interface-trap energy level, no
spread in time constant - gives a Lorentzian peak [[3]](#ref3):

$$
\frac{G_p}{\omega} = qD_{it}A\,\frac{\omega\tau_{it}}{1+(\omega\tau_{it})^2}
\tag{17}
$$

peaking exactly at $\omega\tau_{it}=1$ with height $qD_{it}A/2$. The time
constant itself is set by the capture rate:

$$
\tau_p = \frac{1}{c_pN_A} \tag{18}
$$

with $c_p=\sigma_pv_{th}$ the capture probability [[3]](#ref3).

**Eq. (17) is a teaching-grade starting point, not the full story.** Real
interface traps are spread over a continuum of energies, and real
surfaces have patch-to-patch band-bending fluctuations - both broaden and
lower the measured peak relative to what Eq. (17) predicts
[[3]](#ref3)[[2]](#ref2). Nicollian and Brews give the full statistical
treatment (a Gaussian average over band-bending fluctuations, their
Eqs. 5.70-5.85); illustrating the effect without implementing that full
double integral, superpose a spread of single-level peaks and watch the
sum broaden and flatten relative to any one of them.

```python
Dit_fit, tau_fit = cgv.fit_dit_from_peak(omega, Gp_over_w, area_cm2)
print(f"D_it (true) = {Dit_true:.3e} cm^-2")
print(f"D_it (fit)  = {Dit_fit:.3e} cm^-2  ({100*abs(Dit_fit-Dit_true)/Dit_true:.1f}% off)")
print(f"tau_it (true) = {tau_it_true:.2e} s")
print(f"tau_it (fit)  = {tau_fit:.2e} s")

# Illustrate broadening: superpose single-level peaks over a spread of tau
omega_grid = 2*np.pi*np.logspace(1, 8, 400)
taus = np.geomspace(tau_it_true/8, tau_it_true*8, 9)
single = cgv.conductance_lorentzian(omega_grid, Dit_true, tau_it_true, area_cm2)
broadened = np.mean(
    [cgv.conductance_lorentzian(omega_grid, Dit_true, t, area_cm2) for t in taus], axis=0)

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogx(omega_grid*tau_it_true, single*1e12, lw=2, label='single level (Eq. 17)')
ax.semilogx(omega_grid*tau_it_true, broadened*1e12, lw=2, ls='--',
            label='spread of time constants\n(illustrates the real, broader peak)')
ax.set_xlabel(r'$\omega\tau_{it}$')
ax.set_ylabel(r'$G_p/\omega$ (pF)')
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()
```

![Output 8](assets/nb/cgv_analysis_08.png)

```text
D_it (true) = 5.000e+11 cm^-2
D_it (fit)  = 5.407e+11 cm^-2  (8.1% off)
tau_it (true) = 2.00e-05 s
tau_it (fit)  = 1.83e-05 s
```

The broadened curve peaks lower and wider than the single-level
model, for the same total trap density - reading $D_{it}$ off Eq. (17)'s
peak-height formula when the real peak is this broadened would
underestimate $D_{it}$. This is exactly the caution Nicollian and Brews
give: the width of a measured $G_p/\omega$ curve is itself information
(about the spread of band-bending fluctuations across the surface, i.e.
about surface uniformity), not just noise to average away.

## 14. What ruins a real measurement, and how do you know?

Three distinct non-idealities. They are commonly confused with each
other, and telling them apart is the actual skill this section teaches.

**Series resistance.** Contact and substrate spreading resistance distort
both $C_m$ and $G_m$, worst in strong accumulation at high frequency,
where the device looks almost purely resistive. Extract it there
[[4]](#ref4):

$$
R_s = \frac{(G/\omega C)^2}{\left[1+(G/\omega C)^2\right]G} \tag{19}
$$

then correct the rest of the sweep:

$$
a_R = G_m - \left(G_m^2+(\omega C_m)^2\right)R_s,\qquad
C_{\text{adj}} = \frac{\left(G_m^2+(\omega C_m)^2\right)C_m}{a_R^2+(\omega C_m)^2}
\tag{20}
$$

Left uncorrected, $R_s$ can make even $C_{\text{ox}}$ itself read low.

**Sweep-rate artefacts.** Already met in Sec. 6 as deep depletion: a bias
swept faster than the minority-carrier generation lifetime never reaches
the equilibrium LF curve.

**Mobile ionic charge.** Sodium or potassium ions drifting under bias (or
worse, under bias-temperature stress) shift $V_{FB}$ between a forward
and a reverse sweep - hysteresis. The shift gives the mobile sheet
density directly:

$$
\Delta N_m = -\frac{C_{\text{ox}}\Delta V_{FB}}{qA} \tag{21}
$$

The minus sign matters: positive mobile ions (Na$^+$, K$^+$) that drift to
the oxide-semiconductor interface between the forward and reverse sweep
make $V_{FB}$ *more negative* on the reverse sweep - the same relationship
as Eq. (9), where a more positive $Q_{\text{eff}}$ pulls $V_{FB}$ down.
Dropping the sign would report a negative ion density for a real,
physically positive contamination level.

**These three do not look alike, and should not be confused.** Series
resistance distorts the curve's shape at a fixed bias history. A
sweep-rate artefact changes the inversion plateau, not the depletion
slope. Mobile-ion hysteresis shifts the *whole curve* sideways, by a
different amount on the way up than on the way down, and is the only one
of the three that depends on sweep direction at all.

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

# --- series resistance ---
Rs_true = 40.0
f_test = 1e6
omega_t = 2*np.pi*f_test
Cm_ideal = C_hf  # reuse the Sec. 6 HF curve as the "ideal" device
Gm_ideal = np.full_like(Cm_ideal, 1e-9)
Zmeas = Rs_true + 1/(Gm_ideal + 1j*omega_t*Cm_ideal)
Ymeas = 1/Zmeas
Gm_meas, Cm_meas = Ymeas.real, Ymeas.imag/omega_t
axes[0].plot(Vg, Cm_ideal*1e12, label='ideal (no $R_s$)')
axes[0].plot(Vg, Cm_meas*1e12, label=f'distorted by $R_s$={Rs_true:.0f} $\Omega$')
axes[0].set_xlabel('$V_G$ (V)'); axes[0].set_ylabel('C (pF)')
axes[0].set_title('series resistance'); axes[0].legend(fontsize=8)

# --- deep depletion (sweep rate) ---
axes[1].plot(Vg, C_hf*1e12, label='equilibrium (HF)')
axes[1].plot(Vg, C_dd*1e12, ls='--', label='swept too fast')
axes[1].set_xlabel('$V_G$ (V)'); axes[1].set_ylabel('C (pF)')
axes[1].set_title('sweep-rate artefact'); axes[1].legend(fontsize=8)

# --- mobile-ion hysteresis ---
delta_VFB = -0.35
C_fwd, C_rev = cgv.synthetic_hysteresis_sweep(Vg, area_cm2, tox_nm, Na, 'p', VFB, delta_VFB)
axes[2].plot(Vg, C_fwd*1e12, label='forward sweep')
axes[2].plot(Vg, C_rev*1e12, label='reverse sweep')
axes[2].set_xlabel('$V_G$ (V)'); axes[2].set_ylabel('C (pF)')
axes[2].set_title('mobile-ion hysteresis'); axes[2].legend(fontsize=8)

plt.tight_layout()
plt.show()

Nm = cgv.mobile_charge_from_hysteresis(delta_VFB, Cox, area_cm2)
print(f"Delta V_FB = {delta_VFB:.3f} V  ->  mobile ion density = {Nm:.3e} cm^-2")
```

![Output 9](assets/nb/cgv_analysis_09.png)

```text
Delta V_FB = -0.350 V  ->  mobile ion density = 7.543e+11 cm^-2
```

## 15. Assumptions and limitations

The models here are deliberately simple. They assume:

- **one dimension and a uniform surface**, with no lateral variation in
  doping, oxide thickness, or trap density under the gate.
- **Boltzmann statistics and complete ionisation**, so the model is not
  reliable at degenerate doping.
- **the depletion approximation**, used only as a cross-check (Eq. 6) -
  the working equations use the exact $Q_{sc}(\phi_s)$ of Eq. (3)
  throughout, but that formula itself still assumes a classical,
  non-degenerate semiconductor.
- **a single-time-constant conductance model** (Eq. 17) rather than the
  full band-bending-fluctuation statistical theory [[3]](#ref3); Sec. 13
  shows what the simplification misses.
- **no quantum-mechanical correction to inversion-layer capacitance**,
  which matters for oxides a few nanometres thick on real modern
  devices - this notebook's 10 nm worked example is comfortably in the
  classical regime, but production CMOS gate oxides are not.
- **small-signal AC amplitude** throughout, so every capacitance and
  conductance is a true differential quantity, not a large-signal
  average.
- **teaching-grade parameter values** (silicon constants, mobility,
  doping, oxide charge) - realistic in shape and scale, not measured on
  any specific device.

They are sufficient to explain what a C-V/G-V sweep contains, why the
same algebra spans a MOS capacitor and a bare Schottky diode, and how a
frequency sweep separates a reactive response from a dissipative one.
Quantitative process-control work needs measured silicon parameters, the
full statistical conductance theory, and - for thin oxides - a
quantum-corrected inversion capacitance model.

**And one thing that is not an assumption but a definition:** every
number in this notebook follows the sign convention set in Sec. 3. Any
comparison with the literature has to start by checking which convention
the other source used - Eq. (12)'s absolute value exists specifically
because this trips people up.

## Summary of equations

| # | Equation | Meaning | Source |
|---|---|---|---|
| (1) | $C_{\text{ox}} = \varepsilon_{\text{ox}}A/t_{\text{ox}}$ | oxide capacitance | [[7]](#ref7) |
| (2) | $\phi_F = V_t\ln(N_A/n_i)$ | bulk Fermi potential | [[7]](#ref7) |
| (3) | $Q_{sc}(\phi_s) = -\varepsilon_sE_s$ | exact semiconductor space charge | [[5]](#ref5)[[9]](#ref9) |
| (4) | $C_s(\phi_s) = -dQ_{sc}/d\phi_s$ | differential semiconductor capacitance | [[5]](#ref5) |
| (5) | $L_D = \sqrt{\varepsilon_sV_t/(qN)}$ | Debye length | [[1]](#ref1) |
| (6) | $W\approx\sqrt{2\varepsilon_s\phi_s/(qN)}$ | depletion approximation (cross-check) | [[7]](#ref7) |
| (7) | $C = C_{\text{ox}}C_s/(C_{\text{ox}}+C_s)$ | series combination | [[2]](#ref2) |
| (8) | $V_G = V_{FB}+\phi_s-Q_{sc}/C_{\text{ox}}$ | charge balance | [[2]](#ref2) |
| (9) | $V_{FB} = \phi_{ms}-Q_{\text{eff}}/C_{\text{ox}}$ | flat-band voltage | [[2]](#ref2) |
| (10) | $Q_{\text{eff}} = C_{\text{ox}}(\phi_{ms}-V_{FB})$ | effective oxide charge | [[2]](#ref2) |
| (11) | $V_T = V_{FB}+s[\ldots]+s\cdot2\lvert\phi_F\rvert$ | threshold voltage | [[6]](#ref6) |
| (12) | $N(W) = 2/(q\varepsilon_sA^2\lvert d(1/C^2)/dV\rvert)$ | doping profile | [[1]](#ref1) |
| (13) | $W = \varepsilon_sA/C$ | depth from measured C | [[1]](#ref1) |
| (14) | $1/C^2 = 2(V_{bi}-V)/(q\varepsilon_sA^2N)$ | Mott-Schottky (no oxide) | [[1]](#ref1) |
| (15) | $Y_m = G_m+j\omega C_m$ | measured admittance | [[3]](#ref3) |
| (16) | $G_p/\omega = \omega C_{\text{ox}}^2G_m/[G_m^2+\omega^2(C_{\text{ox}}-C_m)^2]$ | series-to-parallel transform | [[3]](#ref3) |
| (17) | $G_p/\omega = qD_{it}A\,\omega\tau_{it}/(1+(\omega\tau_{it})^2)$ | single-level conductance peak | [[3]](#ref3) |
| (18) | $\tau_p = 1/(c_pN_A)$ | capture time constant | [[3]](#ref3) |
| (19) | $R_s = (G/\omega C)^2/\{[1+(G/\omega C)^2]G\}$ | series resistance (from accumulation) | [[4]](#ref4) |
| (20) | $C_{\text{adj}}, G_{\text{adj}}$ from $a_R$ | series-resistance correction | [[4]](#ref4) |
| (21) | $\Delta N_m = -C_{\text{ox}}\Delta V_{FB}/(qA)$ | mobile ion density from hysteresis | [[2]](#ref2) |

<a id="references"></a>

## References

<a id="ref1"></a>
**[1]** P. A. Barnes, *Capacitance-Voltage (C-V) Characterization of
Semiconductors*, in Characterization of Materials, Vol. 1 (Wiley, 2012).
doi:[10.1002/0471266965.com038](https://doi.org/10.1002/0471266965.com038)

<a id="ref2"></a>
**[2]** E. H. Nicollian and J. R. Brews, *MOS (Metal Oxide Semiconductor)
Physics and Technology* (Wiley, New York, 1982). ISBN 978-0471082227.

<a id="ref3"></a>
**[3]** E. H. Nicollian and A. Goetzberger, *The Si-SiO2 Interface -
Electrical Properties as Determined by the Metal-Insulator-Silicon
Conductance Technique*, Bell Syst. Tech. J. **46**, 1055-1133 (1967).
doi:[10.1002/j.1538-7305.1967.tb01727.x](https://doi.org/10.1002/j.1538-7305.1967.tb01727.x)

<a id="ref4"></a>
**[4]** J. D. Wiley and G. L. Miller, *Series resistance effects in
semiconductor CV profiling*, IEEE Trans. Electron Devices **22**, 265-272
(1975).
doi:[10.1109/T-ED.1975.18109](https://doi.org/10.1109/T-ED.1975.18109)

<a id="ref5"></a>
**[5]** A. S. Grove, E. H. Snow, B. E. Deal and C. T. Sah, *Simple
physical model for the space-charge capacitance of
metal-oxide-semiconductor structures*, J. Appl. Phys. **35**, 2458-2460
(1964).
doi:[10.1063/1.1713760](https://doi.org/10.1063/1.1713760)

<a id="ref6"></a>
**[6]** L. M. Terman, *An investigation of surface states at a
silicon/silicon oxide interface employing metal-oxide-silicon diodes*,
Solid-State Electronics **5**, 285-299 (1962).
doi:[10.1016/0038-1101(62)90111-9](https://doi.org/10.1016/0038-1101(62)90111-9)

<a id="ref7"></a>
**[7]** S. M. Sze and K. K. Ng, *Physics of Semiconductor Devices*, 3rd
ed. (Wiley, Hoboken, 2007).
doi:[10.1002/0470068329](https://doi.org/10.1002/0470068329)

<a id="ref8"></a>
**[8]** D. K. Schroder, *Semiconductor Material and Device
Characterization*, 3rd ed. (Wiley, Hoboken, 2006).
doi:[10.1002/0471749095](https://doi.org/10.1002/0471749095)

<a id="ref9"></a>
**[9]** C. G. B. Garrett and W. H. Brattain, *Physical theory of
semiconductor surfaces*, Phys. Rev. **99**, 376-387 (1955).
doi:[10.1103/PhysRev.99.376](https://doi.org/10.1103/PhysRev.99.376)

<a id="ref10"></a>
**[10]** C. N. Berglund, *Surface states at steam-grown silicon-silicon
dioxide interfaces*, IEEE Trans. Electron Devices **13**, 701-705 (1966).
doi:[10.1109/T-ED.1966.15827](https://doi.org/10.1109/T-ED.1966.15827)
