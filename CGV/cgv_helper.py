"""
cgv_helper.py
-------------
Helper functions for the capacitance-voltage / conductance-voltage (C-V /
G-V) teaching notebook, `cgv_analysis.ipynb`. Equation numbers quoted in
the docstrings below, e.g. "Eq. (12)", refer to the numbered equations in
that notebook.

Scope
-----
These functions implement a 1-D electrostatic description of a MOS
capacitor (metal gate / oxide / doped silicon substrate) under a DC bias
with a small AC probe signal, following exactly the physics used to
extract the parameters a real C-V/G-V measurement reports: oxide
capacitance, flat-band and threshold voltage, a doping profile, and an
interface-trap density from the conductance method. A Schottky (no-oxide)
diode is treated as the same electrostatics with C_ox removed.

They are written for teaching, not for fab-line metrology:
  - the surface is uniform and one-dimensional,
  - Boltzmann statistics throughout; dopants are fully ionised,
  - the "deep-depletion" and interface-trap models are simplified,
    single-mechanism pictures, not the full transient/statistical theory,
  - quantum-mechanical corrections to inversion-layer capacitance
    (important for oxides a few nm thick on real modern devices) are
    ignored.

For quantitative work, use tabulated silicon parameters, a full transient
minority-carrier generation model, and the full band-bending-fluctuation
statistical theory of the conductance method (Nicollian & Brews, Ch. 5).

UNITS - read this before using anything below
------------------------------------------------
Capacitance and charge quantities that are physically "per unit area" -
C_ox, C_s, C (the total measured C-V capacitance), and Q_sc, Q_eff -
are ALL handled in this module AS PER-UNIT-AREA quantities: F/cm^2 and
C/cm^2. Every C-V / G-V curve function therefore takes an explicit
`area_cm2` and returns TOTAL capacitance in F (the pF-nF scale a real
meter reports), by multiplying the per-area physics by area at the very
end. If you need the per-area number directly, divide the returned curve
by `area_cm2`.

Energies and work functions are in eV. Potentials/voltages are in V.
Doping and carrier densities are in cm^-3. D_it is in cm^-2 eV^-1. Oxide
thickness is in nm. Frequencies are given as angular frequency omega
(rad/s) to fit units directly into the admittance formulas; pass
omega = 2*pi*f.

SIGN CONVENTION - read this before using anything below
-----------------------------------------------------------
This module treats a p-type substrate as the default worked example
(matching the Keithley/Tektronix 4200A-SCS application note). Everything
here also works for n-type via `dopant_type`.

  * The surface potential phi_s is the band bending at the
    semiconductor-oxide interface, relative to the neutral bulk, with
    phi_s = 0 at flat band.
  * For p-type: phi_s < 0 is ACCUMULATION, phi_s > 0 is DEPLETION, and
    phi_s > 2*phi_F is (strong) INVERSION, where
    phi_F = V_t*ln(N_A/n_i) > 0 is the bulk Fermi potential (Eq. 2). For
    n-type, phi_F < 0 in this convention and every regime sign flips.
  * Increasing gate bias V_G moves a p-type surface from accumulation
    towards inversion; for n-type it is the reverse.
  * The space-charge density Q_sc(phi_s) (Eq. 3) is SIGNED: positive
    when the surface carries net positive charge (accumulated holes on
    p-type), negative for net negative charge (ionised acceptors in
    depletion on p-type).

References (see the notebook's References section for full citations)
--------------------------------------------------------------------------
  Barnes, P. A., "Capacitance-Voltage (C-V) Characterization of
    Semiconductors", in Characterization of Materials (Wiley, 2012).
    doi:10.1002/0471266965.com038                                [Barnes]
  Nicollian, E. H. and Brews, J. R., MOS (Metal Oxide Semiconductor)
    Physics and Technology (Wiley, 1982).                           [NB]
  Nicollian, E. H. and Goetzberger, A., "The Si-SiO2 Interface -
    Electrical Properties as Determined by the Metal-Insulator-Silicon
    Conductance Technique", Bell Syst. Tech. J. 46, 1055 (1967).   [NG67]
  Wiley, J. D. and Miller, G. L., "Series resistance effects in
    semiconductor C-V profiling", IEEE Trans. Electron Dev. ED-22, 265
    (1975).                                                       [WM75]
  Grove, A. S., Snow, E. H., Deal, B. E. and Sah, C. T., "Simple
    physical model for the space-charge capacitance of
    metal-oxide-semiconductor structures", J. Appl. Phys. 35, 2458
    (1964).                                                      [GSDS64]
"""

import numpy as np
from scipy.optimize import brentq, curve_fit, minimize_scalar
from scipy.signal import savgol_filter

_trapz = getattr(np, "trapezoid", None) or np.trapz

# ---------------------------------------------------------------------------
# Physical constants (matches KPSPV/kpspv_helper.py and TLM/tlm_helper.py
# so numbers are directly comparable across notebooks)
# ---------------------------------------------------------------------------
Q = 1.602176634e-19       # elementary charge, C
KB = 1.380649e-23         # Boltzmann constant, J/K
EPS0 = 8.8541878128e-14   # vacuum permittivity, F/cm

T_ROOM = 300.0            # K

EG_SI_EV = 1.12           # Si band gap, eV
CHI_SI_EV = 4.05          # Si electron affinity, eV
K_SI = 11.7               # Si relative permittivity
EPS_SI = K_SI * EPS0      # F/cm
NI_SI = 9.65e9            # Si intrinsic carrier density, cm^-3
VTH_SI = 1.1e7            # thermal velocity, cm/s

K_SIO2 = 3.9
K_SIN = 7.5
K_AL2O3 = 9.0

WORK_FUNCTION_EV = {
    'Al': 4.10,
    'Au': 5.10,
    'poly-Si (n+)': 4.05,
    'poly-Si (p+)': 5.17,
    'W': 4.55,
}


def thermal_voltage(T=T_ROOM):
    """Thermal voltage V_t = kT/q, volts."""
    return KB * T / Q


def _sign_dope(dopant_type):
    """+1 for p-type, -1 for n-type - the convention used throughout."""
    if dopant_type == 'p':
        return 1.0
    if dopant_type == 'n':
        return -1.0
    raise ValueError("dopant_type must be 'p' or 'n'")


# ===========================================================================
# PART 1 - the MOS structure: oxide capacitance, work functions
# ===========================================================================
def oxide_capacitance_per_area(tox_nm, k_ox=K_SIO2):
    """Oxide capacitance per unit area, C_ox = eps_ox/t_ox - Eq. (1),
    F/cm^2."""
    return (k_ox * EPS0) / (tox_nm * 1e-7)


def oxide_capacitance(area_cm2, tox_nm, k_ox=K_SIO2):
    """Total oxide capacitance C_ox = eps_ox*A/t_ox - Eq. (1), F."""
    return oxide_capacitance_per_area(tox_nm, k_ox) * area_cm2


def fermi_potential(doping_cm3, dopant_type, T=T_ROOM):
    """
    Bulk Fermi potential phi_F = V_t*ln(N/n_i) - Eq. (2): positive for
    p-type, negative for n-type. Anchors the inversion condition
    phi_s = 2*phi_F and the threshold-voltage formula, Eq. (11).
    """
    Vt = thermal_voltage(T)
    return _sign_dope(dopant_type) * Vt * np.log(doping_cm3 / NI_SI)


def debye_length_cm(doping_cm3, T=T_ROOM):
    """Extrinsic Debye length L_D = sqrt(eps_s*V_t/(q*N)) - Eq. (5). Sets
    the depth resolution of C-V doping profiling and the flat-band
    capacitance [Barnes Eq. 27]."""
    Vt = thermal_voltage(T)
    return np.sqrt(EPS_SI * Vt / (Q * doping_cm3))


def work_function_ms(phi_metal_eV, chi_s_eV, Eg_eV, phi_F_eV):
    """Metal-semiconductor work function difference
    phi_ms = phi_metal - (chi_s + Eg/2 - phi_F), after [NB]. phi_F_eV is
    SIGNED per this module's convention."""
    return phi_metal_eV - (chi_s_eV + Eg_eV / 2.0 - phi_F_eV)


# ===========================================================================
# PART 2 - the semiconductor surface: exact space charge and its
# differential capacitance, valid through accumulation/depletion/inversion
# ===========================================================================
def surface_field(phi_s, doping_cm3, dopant_type, T=T_ROOM):
    """
    Exact semiconductor surface field (not the depletion approximation),
    valid in accumulation, depletion AND inversion [GSDS64]:

        E_s = 2*sign(phi_s)*sqrt(q*n_i*V_t/eps_s)
              * sqrt( cosh((phi_s-phi_F)/V_t) + (phi_s/V_t)*sinh(phi_F/V_t)
                      - cosh(phi_F/V_t) )

    V/cm. phi_s = 0 (flat band) gives E_s = 0 exactly.
    """
    phi_s = np.asarray(phi_s, dtype=float)
    Vt = thermal_voltage(T)
    phiF = fermi_potential(doping_cm3, dopant_type, T)
    arg = (np.cosh((phi_s - phiF) / Vt)
           + (phi_s / Vt) * np.sinh(phiF / Vt)
           - np.cosh(phiF / Vt))
    arg = np.clip(arg, 0.0, None)
    pref = 2.0 * np.sqrt(Q * NI_SI * Vt / EPS_SI)
    return np.sign(phi_s) * pref * np.sqrt(arg)


def space_charge_density(phi_s, doping_cm3, dopant_type, T=T_ROOM):
    """Semiconductor space-charge density Q_sc(phi_s), SIGNED, C/cm^2 -
    Eq. (3). Q_sc = -eps_s*E_s (Gauss's law)."""
    return -EPS_SI * surface_field(phi_s, doping_cm3, dopant_type, T)


def semiconductor_capacitance(phi_s, doping_cm3, dopant_type, T=T_ROOM):
    """
    Differential semiconductor capacitance per unit area,
    C_s(phi_s) = -dQ_sc/dphi_s - Eq. (4), F/cm^2. Central difference of
    the exact Q_sc(phi_s) of Eq. (3): a DIFFERENTIAL quantity, not DC
    charge/DC potential.

    Q_sc(phi_s) is a MONOTONICALLY DECREASING function of phi_s (more
    band bending always drives the semiconductor charge more negative,
    e.g. p-type: positive in accumulation, negative in depletion), so the
    physical (positive) capacitance is the NEGATIVE of dQ_sc/dphi_s - the
    gate charge -Q_sc is what increases with phi_s.
    """
    phi_s = np.asarray(phi_s, dtype=float)
    h = 1e-4  # V
    q1 = space_charge_density(phi_s + h, doping_cm3, dopant_type, T)
    q2 = space_charge_density(phi_s - h, doping_cm3, dopant_type, T)
    return -(q1 - q2) / (2 * h)


def depletion_width_depletion_approx_cm(phi_s, doping_cm3):
    """Depletion-approximation width W = sqrt(2*eps_s*|phi_s|/(q*N)) -
    Eq. (6). Analytic cross-check on Eq. (3)/(4); valid only in
    depletion, 0 < phi_s < ~2*phi_F."""
    return np.sqrt(2 * EPS_SI * np.abs(phi_s) / (Q * doping_cm3))


# ===========================================================================
# PART 3 - the ideal C-V curve (all curve functions return TOTAL
# capacitance in F, i.e. per-area physics x area_cm2)
# ===========================================================================
def total_capacitance_series_per_area(Cox_pa, Cs_pa):
    """Series combination, per unit area: C = Cox*Cs/(Cox+Cs) - Eq. (7).
    F/cm^2."""
    return Cox_pa * Cs_pa / (Cox_pa + Cs_pa)


def surface_potential_from_bias(Vg, Cox_pa, VFB, doping_cm3, dopant_type,
                                 T=T_ROOM, phi_s_max_factor=6.0):
    """
    Solve the charge-balance relation for phi_s at each gate bias -
    Eq. (8):  V_G = V_FB + phi_s - Q_sc(phi_s)/Cox_pa   (per unit area
    throughout, so area cancels). Returns phi_s (V), same shape as Vg.
    """
    scalar_input = np.isscalar(Vg) or np.asarray(Vg).ndim == 0
    Vg = np.atleast_1d(np.asarray(Vg, dtype=float))
    phiF = fermi_potential(doping_cm3, dopant_type, T)
    hi = phi_s_max_factor * abs(phiF) + 1.0
    out = np.empty_like(Vg)

    def resid(phi_s, vg):
        Qsc = space_charge_density(phi_s, doping_cm3, dopant_type, T)
        return (VFB + phi_s - Qsc / Cox_pa) - vg

    for i, vg in enumerate(Vg):
        out[i] = brentq(resid, -hi, hi, args=(vg,), xtol=1e-12, rtol=1e-13)
    return out[0] if scalar_input else out


def flatband_voltage(phi_ms_eV, Qeff_C_cm2, Cox_pa):
    """Flat-band voltage V_FB = phi_ms - Q_eff/C_ox - Eq. (9), the MOS
    form of [NB] Eq. 10.10. Cox_pa is per unit area (F/cm^2)."""
    return phi_ms_eV - Qeff_C_cm2 / Cox_pa


def effective_oxide_charge(phi_ms_eV, VFB, Cox_pa):
    """Effective oxide charge Q_eff = Cox*(phi_ms - V_FB) - Eq. (10),
    rearranging Eq. (9). Returns C/cm^2."""
    return Cox_pa * (phi_ms_eV - VFB)


def flatband_capacitance_per_area(Cox_pa, doping_cm3, T=T_ROOM):
    """Flat-band capacitance C_FB = series(Cox, eps_s/L_D) - used to find
    V_FB on a measured curve (Tektronix 4200A-SCS application note).
    F/cm^2."""
    Cs_fb = EPS_SI / debye_length_cm(doping_cm3, T)
    return total_capacitance_series_per_area(Cox_pa, Cs_fb)


def threshold_voltage(Cox_pa, doping_cm3, dopant_type, VFB, T=T_ROOM):
    """
    Threshold voltage - Eq. (11):
        V_T = V_FB + s*sqrt(4*eps_s*q*N*|phi_F|)/Cox + s*2*|phi_F|
    s = +1 (p-type, V_T > V_FB) or -1 (n-type, V_T < V_FB), after the
    Tektronix 4200A-SCS application note.
    """
    s = _sign_dope(dopant_type)
    phiF = fermi_potential(doping_cm3, dopant_type, T)
    depletion_term = np.sqrt(4 * EPS_SI * Q * doping_cm3 * abs(phiF)) / Cox_pa
    return VFB + s * depletion_term + s * 2 * abs(phiF)


def _phi_s_min_capacitance(doping_cm3, dopant_type, T=T_ROOM):
    """
    The surface potential at which the EXACT equilibrium C_s(phi_s) of
    Eq. (4) is smallest. This is close to, but not exactly at, the
    textbook "onset of strong inversion" phi_s = 2*phi_F: the exact
    formula's minority-carrier term only starts winning against the
    still-shrinking depletion term a little past that point. Used to
    clamp the HF curve at its true minimum instead of at the 2*phi_F
    approximation, which would otherwise leave a small non-physical step
    in the HF curve where clamping begins.
    """
    s = _sign_dope(dopant_type)
    phiF = fermi_potential(doping_cm3, dopant_type, T)
    lo, hi = sorted([0.2 * s * phiF, 4.0 * s * phiF])

    def log_Cs(phi_s):
        return np.log(semiconductor_capacitance(np.array([phi_s]),
                                                  doping_cm3, dopant_type, T)[0])

    res = minimize_scalar(log_Cs, bounds=(lo, hi), method="bounded",
                           options={"xatol": 1e-6})
    return res.x


def hf_cv_curve(Vg, area_cm2, tox_nm, doping_cm3, dopant_type, VFB,
                T=T_ROOM, k_ox=K_SIO2):
    """
    High-frequency C-V curve - Eqs. (6)-(7). Minority carriers cannot
    follow the AC signal once the surface reaches strong inversion, so
    C_s is frozen at its value at the true minimum of the equilibrium
    C_s(phi_s) curve (see `_phi_s_min_capacitance`; close to, but not
    exactly, the textbook phi_s = 2*phi_F) for any more strongly-inverting
    phi_s. Returns TOTAL capacitance, F.
    """
    Cox_pa = oxide_capacitance_per_area(tox_nm, k_ox)
    phi_s_min = _phi_s_min_capacitance(doping_cm3, dopant_type, T)
    s = _sign_dope(dopant_type)
    phi_s = surface_potential_from_bias(Vg, Cox_pa, VFB, doping_cm3,
                                         dopant_type, T)
    phi_s_hf = np.where(s * phi_s > s * phi_s_min, phi_s_min, phi_s)
    Cs_pa = semiconductor_capacitance(phi_s_hf, doping_cm3, dopant_type, T)
    return total_capacitance_series_per_area(Cox_pa, Cs_pa) * area_cm2


def lf_cv_curve(Vg, area_cm2, tox_nm, doping_cm3, dopant_type, VFB,
                 T=T_ROOM, k_ox=K_SIO2):
    """
    Low-frequency / quasi-static C-V curve - Eqs. (6)-(7). Minority
    carriers are assumed to fully equilibrate at every bias step: the
    exact Q_sc(phi_s) of Eq. (3) already contains the equilibrium
    inversion charge, so no clamping is applied. Returns TOTAL
    capacitance, F.
    """
    Cox_pa = oxide_capacitance_per_area(tox_nm, k_ox)
    phi_s = surface_potential_from_bias(Vg, Cox_pa, VFB, doping_cm3,
                                         dopant_type, T)
    Cs_pa = semiconductor_capacitance(phi_s, doping_cm3, dopant_type, T)
    return total_capacitance_series_per_area(Cox_pa, Cs_pa) * area_cm2


def deep_depletion_cv(Vg, area_cm2, tox_nm, doping_cm3, dopant_type, VFB,
                       T=T_ROOM, k_ox=K_SIO2):
    """
    Deep-depletion C-V curve - the sweep-too-fast limit: minority
    carriers get NO chance to generate, so depletion keeps growing past
    its equilibrium turning point at phi_s = 2*phi_F instead of
    saturating there. Modelled with the depletion-approximation charge
    -sign(phi_s)*sqrt(2*eps_s*q*N*|phi_s|) continued for phi_s beyond
    2*phi_F (majority-carrier depletion charge only, no inversion
    charge), stitched to the exact Eq. (3) result below that point so the
    accumulation/depletion branch is unchanged. Returns TOTAL
    capacitance, F.
    """
    Cox_pa = oxide_capacitance_per_area(tox_nm, k_ox)
    phi_s_inv = _phi_s_min_capacitance(doping_cm3, dopant_type, T)
    phiF = fermi_potential(doping_cm3, dopant_type, T)
    s = _sign_dope(dopant_type)

    def Qsc_deep(phi_s):
        phi_s = np.atleast_1d(np.asarray(phi_s, dtype=float))
        Q_exact = space_charge_density(phi_s, doping_cm3, dopant_type, T)
        beyond = s * phi_s > s * phi_s_inv
        Q_dep_only = -np.sign(phi_s) * np.sqrt(
            2 * EPS_SI * Q * doping_cm3 * np.abs(phi_s))
        return np.where(beyond, Q_dep_only, Q_exact)

    Vg = np.atleast_1d(np.asarray(Vg, dtype=float))
    hi = 6 * abs(phiF) + 1.0
    phi_s = np.empty_like(Vg)
    for i, vg in enumerate(Vg):
        resid = lambda p: (VFB + p - Qsc_deep(p)[0] / Cox_pa) - vg
        phi_s[i] = brentq(resid, -hi, hi, xtol=1e-12, rtol=1e-13)

    h = 1e-4
    Cs_pa = (Qsc_deep(phi_s + h) - Qsc_deep(phi_s - h)) / (2 * h)
    return total_capacitance_series_per_area(Cox_pa, Cs_pa) * area_cm2


# ===========================================================================
# PART 4 - doping profile from C(V), and the Mott-Schottky (no-oxide) limit
# ===========================================================================
def doping_profile_from_cv(V, C, area_cm2, eps_s=EPS_SI, smooth=False):
    """
    Doping density vs. depth from a measured HF C-V curve - Eqs. (11)-(12)
    [Barnes Eqs. 20-21; Tektronix 4200A-SCS application note]:

        N(W) = -2 / (q*eps_s*A^2 * d(1/C^2)/dV),   W = eps_s*A/C

    C is TOTAL capacitance (F), V the gate/reverse bias (V), same length.
    If smooth=True, 1/C^2 is smoothed with a Savitzky-Golay filter before
    differentiating (a regularised derivative, contrasted in Sec. 9 with
    the raw finite difference) - both are legitimate; the raw one is
    deliberately noisier, which is the point.
    """
    V = np.asarray(V, dtype=float)
    C = np.asarray(C, dtype=float)
    order = np.argsort(V)
    V, C = V[order], C[order]
    inv_c2 = 1.0 / C**2
    if smooth and len(V) >= 7:
        window = min(11, len(V) - (1 - len(V) % 2))
        if window % 2 == 0:
            window -= 1
        window = max(window, 5)
        inv_c2 = savgol_filter(inv_c2, window, 2)
    d_inv_c2_dV = np.gradient(inv_c2, V)
    # |d(1/C^2)/dV| rather than a signed derivative: the textbook formula
    # is usually written for a bias convention where increasing voltage
    # means increasing REVERSE bias (so 1/C^2 falls as V rises); the
    # sweep convention used elsewhere in this module has increasing V_G
    # deepen depletion directly for p-type (1/C^2 RISES as V rises), the
    # opposite sign. N(W) itself is never negative, so taking the
    # magnitude of the slope sidesteps having to track that sign
    # convention here - see Sec. 9's discussion of this exact trap.
    N = 2.0 / (Q * eps_s * area_cm2**2 * np.abs(d_inv_c2_dV))
    W = eps_s * area_cm2 / C
    return W, N


def mott_schottky_fit(V, C, area_cm2, eps_s=EPS_SI):
    """
    Mott-Schottky analysis - Eq. (14) [Barnes Eq. 14]:

        1/C^2 = 2*(V_bi - V) / (q*eps_s*A^2*N)

    linear in V. Returns (N, N_err, Vbi, Vbi_err) from a weighted linear
    fit of 1/C^2 vs V (slope -> N, intercept/slope -> Vbi). This is the
    C_ox -> infinity limit of Eqs. (6),(11) - see Sec. 10.
    """
    V = np.asarray(V, dtype=float)
    C = np.asarray(C, dtype=float)
    y = 1.0 / C**2
    A = np.vstack([V, np.ones_like(V)]).T
    (slope, intercept), residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    n = len(V)
    y_fit = A @ np.array([slope, intercept])
    dof = max(n - 2, 1)
    resid_var = np.sum((y - y_fit) ** 2) / dof
    cov = resid_var * np.linalg.inv(A.T @ A)
    slope_err = np.sqrt(cov[0, 0])
    intercept_err = np.sqrt(cov[1, 1])

    N = -2.0 / (Q * eps_s * area_cm2**2 * slope)
    N_err = abs(N) * abs(slope_err / slope)
    Vbi = -intercept / slope
    Vbi_err = abs(Vbi) * np.sqrt((intercept_err / intercept) ** 2
                                  + (slope_err / slope) ** 2)
    return N, N_err, Vbi, Vbi_err


# ===========================================================================
# PART 5 - the conductance method (interface trap density)
# ===========================================================================
def admittance_to_parallel(Cm, Gm, omega, Cox_total):
    """
    Convert a MEASURED admittance (Cm, Gm, at angular frequency omega,
    through the series oxide capacitance Cox_total) to the equivalent
    parallel Gp/omega of the semiconductor branch alone - Eq. (16)
    [NB Eq. 5.77]:

        Gp/omega = omega*Cox^2*Gm / (Gm^2 + omega^2*(Cox-Cm)^2)

    All capacitances/conductances here are TOTAL (F, S), matching what an
    LCR meter reports; Cox_total must be measured in strong accumulation.
    """
    return (omega * Cox_total**2 * Gm) / (Gm**2 + (omega * (Cox_total - Cm))**2)


def conductance_lorentzian(omega, Dit, tau_it, area_cm2):
    """
    Single-time-constant (no band-bending dispersion) conductance model -
    Eq. (17), the simplest limit of the full statistical theory
    [NB Ch. 5, Fig. 5.14 curve (c)]:

        Gp/omega = q*Dit*A * omega*tau_it / (1 + (omega*tau_it)^2)

    Peaks at omega*tau_it = 1 with peak height q*Dit*A/2. The REAL
    Gp/omega peak is lower and broader than this because of the
    band-bending fluctuations this simple model omits (Sec. 13; see
    NB Eqs. 5.70-5.85 for the full treatment).
    """
    x = omega * tau_it
    return Q * Dit * area_cm2 * x / (1 + x**2)


def fit_dit_from_peak(omega_array, Gp_over_omega_array, area_cm2):
    """
    Extract (Dit, tau_it) from a measured Gp/omega vs omega sweep at
    fixed bias, by a nonlinear least-squares fit of the single-level
    Lorentzian (Eq. 17) - the simple limit of [NB Eq. 5.81] (sigma_s -> 0)
    - to every point in the sweep. Uses the peak condition omega*tau_it=1,
    (Gp/omega)_max=q*Dit*A/2 only to seed the fit's initial guess; reading
    Dit off a single sample point (the literal peak-height formula) is not
    used for the returned value because it is far more sensitive to
    per-point measurement noise than a fit using the whole curve is.
    Returns (Dit, tau_it).
    """
    omega_array = np.asarray(omega_array, dtype=float)
    Gp_over_omega_array = np.asarray(Gp_over_omega_array, dtype=float)

    # A single noisy sample right at the peak is a bad seed - it can pull
    # curve_fit into a nearby local minimum instead of the true one (the
    # Lorentzian is sharp on a log-omega grid). A light Savitzky-Golay
    # smoothing pass only sets the initial guess; the fit itself still
    # runs against the raw (unsmoothed) data.
    window = min(9, len(Gp_over_omega_array))
    if window % 2 == 0:
        window -= 1
    if window >= 5:
        smoothed = savgol_filter(Gp_over_omega_array, window_length=window,
                                  polyorder=2)
    else:
        smoothed = Gp_over_omega_array
    i_peak = np.argmax(smoothed)
    tau0 = 1.0 / omega_array[i_peak]
    Dit0 = 2 * smoothed[i_peak] / (Q * area_cm2)

    def _model(omega, Dit, tau_it):
        x = omega * tau_it
        return Q * Dit * area_cm2 * x / (1 + x**2)

    popt, _ = curve_fit(_model, omega_array, Gp_over_omega_array,
                         p0=[Dit0, tau0], bounds=([0.0, 0.0], [np.inf, np.inf]))
    return popt[0], popt[1]


def time_constant_from_capture(sigma_cm2, doping_cm3, v_th_cm_s=VTH_SI):
    """Capture time constant tau_p = 1/(c_p*N_A) [NB Eq. 5.74b], with
    c_p = sigma*v_th the capture probability. Seconds."""
    c_p = sigma_cm2 * v_th_cm_s
    return 1.0 / (c_p * doping_cm3)


# ===========================================================================
# PART 6 - non-idealities: series resistance and mobile charge
# ===========================================================================
def extract_series_resistance(Cm_acc, Gm_acc, omega):
    """
    Series resistance from a measurement taken deep in strong
    accumulation, where the device looks almost purely resistive - Eq. (19) (Keithley/Tektronix 4200A-SCS application note, after [WM75]):

        Rs = (G/(omega*C))^2 / { [1 + (G/(omega*C))^2] * G }

    Ohms.
    """
    x = Gm_acc / (omega * Cm_acc)
    return x**2 / ((1 + x**2) * Gm_acc)


def series_resistance_correction(Cm, Gm, omega, Rs):
    """
    Series-resistance-corrected parallel capacitance and conductance -
    Eq. (20) (Keithley/Tektronix 4200A-SCS application note, after
    [WM75]):

        a_R = Gm - (Gm^2 + (omega*Cm)^2) * Rs
        C_adj = (Gm^2 + (omega*Cm)^2)*Cm / (a_R^2 + (omega*Cm)^2)
        G_adj = (Gm^2 + (omega*Cm)^2)*a_R / (a_R^2 + (omega*Cm)^2)

    Returns (C_adj, G_adj).
    """
    a_R = Gm - (Gm**2 + (omega * Cm) ** 2) * Rs
    denom = a_R**2 + (omega * Cm) ** 2
    C_adj = (Gm**2 + (omega * Cm) ** 2) * Cm / denom
    G_adj = (Gm**2 + (omega * Cm) ** 2) * a_R / denom
    return C_adj, G_adj


def mobile_charge_from_hysteresis(delta_VFB, Cox_total, area_cm2):
    """Mobile ion sheet density from a forward/reverse V_FB hysteresis -
    Eq. (21): Delta N_m = -Cox*Delta V_FB/(q*A). Positive mobile ions
    (e.g. Na+) that migrate to the oxide-semiconductor interface between
    the forward and reverse sweep make V_FB more negative on the reverse
    sweep (same sign relationship as Eq. 9/10: V_FB = phi_ms - Q_eff/Cox),
    so a negative delta_VFB (reverse - forward) must give a positive ion
    density - the minus sign is required, not optional. Returns cm^-2."""
    return -Cox_total * delta_VFB / (Q * area_cm2)


# ===========================================================================
# PART 7 - synthetic measurements with realistic, deliberate defects
# ===========================================================================
def synthetic_hf_cv(Vg, area_cm2, tox_nm, doping_cm3, dopant_type, VFB,
                     T=T_ROOM, k_ox=K_SIO2, noise_frac=0.0, series_R=0.0,
                     omega=2 * np.pi * 1e6, seed=0):
    """HF C-V with optional multiplicative noise and series-resistance
    distortion (Eqs. 6-7, 19; deliberate defects for Sec. 9/14)."""
    C_true = hf_cv_curve(Vg, area_cm2, tox_nm, doping_cm3, dopant_type,
                          VFB, T, k_ox)
    G_true = np.zeros_like(C_true)  # ideal capacitor: no loss before Rs
    if series_R:
        # Realistic MOS device: assign a small intrinsic Gp so the
        # circuit is non-degenerate, then distort with Rs.
        G_true = 1e-9 * np.ones_like(C_true)
        denom = (1 + series_R * G_true) ** 2 + (series_R * omega * C_true) ** 2
        C_meas = C_true / denom
        G_meas = (G_true + series_R * omega**2 * C_true**2) / denom
    else:
        C_meas = C_true
    rng = np.random.default_rng(seed)
    if noise_frac:
        C_meas = C_meas * (1 + rng.normal(0, noise_frac, size=C_meas.shape))
    return C_meas


def synthetic_mott_schottky(V, N, Vbi, area_cm2, eps_s=EPS_SI,
                             noise_frac=0.0, seed=0):
    """Synthetic 1/C^2-linear Schottky data - Eq. (14), with optional
    noise (Sec. 10)."""
    inv_c2 = 2 * (Vbi - V) / (Q * eps_s * area_cm2**2 * N)
    rng = np.random.default_rng(seed)
    if noise_frac:
        inv_c2 = inv_c2 * (1 + rng.normal(0, noise_frac, size=inv_c2.shape))
    return 1.0 / np.sqrt(np.clip(inv_c2, 1e-30, None))


def synthetic_conductance_sweep(freq_array, Dit, tau_it, area_cm2, tox_nm,
                                 k_ox=K_SIO2, doping_cm3=1e16,
                                 dopant_type='p', VFB=-0.9, Vg_bias=0.2,
                                 T=T_ROOM, noise_frac=0.0, seed=0):
    """
    Synthetic measured (Cm, Gm) at a fixed bias across a frequency sweep,
    for a MOS capacitor whose interface traps follow the single-level
    model, Eq. (17) - i.e. what a real conductance-method measurement
    looks like BEFORE the Eq. (16) transform is applied. Returns
    (Cm, Gm) arrays, same length as freq_array (Hz).
    """
    omega = 2 * np.pi * np.asarray(freq_array, dtype=float)
    Cox_pa = oxide_capacitance_per_area(tox_nm, k_ox)
    phi_s = surface_potential_from_bias(Vg_bias, Cox_pa, VFB, doping_cm3,
                                         dopant_type, T)
    Cs_pa = semiconductor_capacitance(phi_s, doping_cm3, dopant_type, T)
    Cs_total = Cs_pa * area_cm2
    Cox_total = Cox_pa * area_cm2
    Gp_over_omega = conductance_lorentzian(omega, Dit, tau_it, area_cm2)
    Gp = Gp_over_omega * omega
    # Semiconductor branch admittance: Ys = jw*Cs + Gp; combine in series
    # with jw*Cox to get the admittance an LCR meter actually reports.
    Ys = Gp + 1j * omega * Cs_total
    Zs = 1.0 / Ys
    Zox = 1.0 / (1j * omega * Cox_total)
    Ym = 1.0 / (Zs + Zox)
    Cm = Ym.imag / omega
    Gm = Ym.real
    rng = np.random.default_rng(seed)
    if noise_frac:
        Cm = Cm * (1 + rng.normal(0, noise_frac, size=Cm.shape))
        Gm = Gm * (1 + rng.normal(0, noise_frac, size=Gm.shape))
    return Cm, Gm


def synthetic_hysteresis_sweep(Vg, area_cm2, tox_nm, doping_cm3,
                                dopant_type, VFB_forward, delta_VFB, T=T_ROOM,
                                k_ox=K_SIO2):
    """Forward- and reverse-swept HF C-V curves separated by a mobile-ion
    V_FB shift - Eq. (21) (Sec. 14). Returns (C_forward, C_reverse)."""
    C_fwd = hf_cv_curve(Vg, area_cm2, tox_nm, doping_cm3, dopant_type,
                         VFB_forward, T, k_ox)
    C_rev = hf_cv_curve(Vg, area_cm2, tox_nm, doping_cm3, dopant_type,
                         VFB_forward + delta_VFB, T, k_ox)
    return C_fwd, C_rev
