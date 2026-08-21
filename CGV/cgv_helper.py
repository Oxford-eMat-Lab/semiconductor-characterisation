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
from functools import lru_cache
from scipy.optimize import brentq, curve_fit, minimize_scalar
from scipy.signal import savgol_filter

_trapz = getattr(np, "trapezoid", None) or np.trapz


@lru_cache(maxsize=16)
def _leggauss_cached(n):
    return np.polynomial.legendre.leggauss(n)

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

    # The depletion-only branch is shifted so that it meets the exact
    # branch exactly at phi_s_inv. Without this the total charge jumps by
    # ~1e-9 C/cm^2 at the stitch, and the numerical derivative below turns
    # that step into a spurious spike in C_s right at threshold.
    _Q_exact_at_stitch = space_charge_density(
        np.atleast_1d(phi_s_inv), doping_cm3, dopant_type, T)[0]
    _Q_dep_at_stitch = -np.sign(phi_s_inv) * np.sqrt(
        2 * EPS_SI * Q * doping_cm3 * abs(phi_s_inv))
    _stitch_offset = _Q_exact_at_stitch - _Q_dep_at_stitch

    def Qsc_deep(phi_s):
        phi_s = np.atleast_1d(np.asarray(phi_s, dtype=float))
        Q_exact = space_charge_density(phi_s, doping_cm3, dopant_type, T)
        beyond = s * phi_s > s * phi_s_inv
        Q_dep_only = -np.sign(phi_s) * np.sqrt(
            2 * EPS_SI * Q * doping_cm3 * np.abs(phi_s)) + _stitch_offset
        return np.where(beyond, Q_dep_only, Q_exact)

    Vg = np.atleast_1d(np.asarray(Vg, dtype=float))
    hi = 6 * abs(phiF) + 1.0
    phi_s = np.empty_like(Vg)
    for i, vg in enumerate(Vg):
        resid = lambda p: (VFB + p - Qsc_deep(p)[0] / Cox_pa) - vg
        phi_s[i] = brentq(resid, -hi, hi, xtol=1e-12, rtol=1e-13)

    # C_s = -dQ_sc/dphi_s. The minus sign is not optional: Q_sc(phi_s) is
    # monotonically decreasing, so the bare derivative is negative, and a
    # negative C_s in the series formula produces a pole wherever
    # C_s = -C_ox - which is exactly the +/-600 pF excursion this function
    # used to draw. Same convention as semiconductor_capacitance().
    h = 1e-4
    Cs_pa = -(Qsc_deep(phi_s + h) - Qsc_deep(phi_s - h)) / (2 * h)
    return total_capacitance_series_per_area(Cox_pa, Cs_pa) * area_cm2


# ===========================================================================
# PART 4 - doping profile from C(V), and the Mott-Schottky (no-oxide) limit
# ===========================================================================
def doping_profile_from_cv(V, C, area_cm2, eps_s=EPS_SI, smooth=False,
                           mask_invalid=True, c_tol=0.02, slope_floor=0.2,
                           cox_F=None):
    """
    Doping density vs. depth from a measured HF C-V curve - Eqs. (11)-(12)
    [Barnes Eqs. 20-21; Tektronix 4200A-SCS application note]:

        N(W) = -2 / (q*eps_s*A^2 * d(1/C^2)/dV)

    with the depth W measured from the SEMICONDUCTOR capacitance:

        W = eps_s*A*(1/C - 1/C_ox)        (MOS, cox_F given)
        W = eps_s*A/C                     (Schottky / no oxide, cox_F=None)

    Pass cox_F (the oxide capacitance of the device, in farads) for a MOS
    capacitor. Leaving it None treats the whole measured capacitance as the
    semiconductor's, which is correct only when there is no insulator in
    series - the Mott-Schottky case of Sec. 10.

    Getting this wrong does not change N, only the depth it is plotted at,
    which is why it survives a uniform-doping test unnoticed: every depth
    is shifted by the constant eps_s*A/C_ox = eps_s*t_ox/eps_ox, about
    30 nm for 10 nm of SiO2 on silicon. On a non-uniform profile that moves
    an implant peak by 30 nm and is a real error.

    C is TOTAL capacitance (F), V the gate/reverse bias (V), same length.
    If smooth=True, 1/C^2 is smoothed with a Savitzky-Golay filter before
    differentiating (a regularised derivative, contrasted in Sec. 9 with
    the raw finite difference) - both are legitimate; the raw one is
    deliberately noisier, which is the point.

    VALIDITY. The formula divides by d(1/C^2)/dV, and that derivative goes
    to zero as soon as the HF curve saturates in inversion: the depletion
    edge is pinned at W_max and stops responding to bias, so there is no
    longer any depth information in the measurement. Dividing by it
    produces an N that diverges - values of 1e18-1e19 cm^-3 on a 3e16
    substrate - which is an artefact, not a doping spike. With
    mask_invalid=True (the default) those points are returned as NaN.

    The test is done on the capacitance rather than on the derivative,
    because that is what a practitioner actually looks at: a point is
    rejected when C has come within c_tol (default 2%) of the minimum
    capacitance reached in the sweep, i.e. when W is within 2% of
    W_max = eps_s*A/C_min. Pass mask_invalid=False to see the raw
    divergence - Sec. 9 does exactly that, deliberately, before showing
    the masked version.
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
    with np.errstate(divide='ignore', invalid='ignore'):
        N = 2.0 / (Q * eps_s * area_cm2**2 * np.abs(d_inv_c2_dV))
        if cox_F is None:
            W = eps_s * area_cm2 / C
        else:
            W = eps_s * area_cm2 * (1.0 / C - 1.0 / cox_F)

    if mask_invalid:
        # two independent symptoms of the same thing - the depletion edge
        # has stopped moving, so there is no depth information left.
        saturated = C <= C.min() * (1.0 + c_tol)
        slope = np.abs(d_inv_c2_dV)
        collapsed = slope < slope_floor * np.median(slope)
        N = np.where(saturated | collapsed, np.nan, N)
        W = np.where(W > 0, W, np.nan)
    return W, N


def cv_from_depth_profile(depth_cm, N_cm3, area_cm2, tox_nm,
                          dopant_type='p', VFB=0.0, k_ox=K_SIO2,
                          eps_s=EPS_SI):
    """
    Forward-model an HF C-V curve from an arbitrary doping profile N(x) -
    the depletion approximation, exactly.

    For a depletion edge at W, integrating Poisson once and twice gives

        Q_dep(W) = q * int_0^W N(x) dx                (C/cm^2)
        phi_s(W) = (q/eps_s) * int_0^W x N(x) dx      (V)
        V_G(W)   = V_FB + s*(phi_s + Q_dep/C_ox)
        C_s(W)   = eps_s / W,   C = series(C_ox, C_s) * A

    with s = +1 for p-type (positive gate bias depletes) and -1 for n-type.
    The second relation is the depth-weighted first moment of the profile,
    which is what makes a *non-uniform* N(x) tractable: phi_s is not simply
    proportional to W^2 any more.

    This replaces the earlier approach of evaluating W = sqrt(2 eps phi/(q N))
    with a *locally* uniform N at each bias. That was wrong in two ways: it
    let W jump discontinuously (by 102% in one bias step across a doping
    step, when a depletion edge must move continuously), and it dropped the
    oxide voltage drop Q_dep/C_ox entirely - 0.13-0.33 V out of a 1.22 V
    sweep. Both fed straight into d(1/C^2)/dV and made the extracted
    profile meaningless.

    Parameters
    ----------
    depth_cm : increasing depth grid, cm. Must start at (or very near) 0.
    N_cm3    : doping at each depth, cm^-3, same length.

    Returns (V_G, C_total) with V_G monotonic by construction.
    """
    x = np.asarray(depth_cm, dtype=float)
    N = np.asarray(N_cm3, dtype=float)
    if x.ndim != 1 or x.shape != N.shape:
        raise ValueError("depth_cm and N_cm3 must be 1-D and the same length")
    if np.any(np.diff(x) <= 0):
        raise ValueError("depth_cm must be strictly increasing")

    Cox_pa = oxide_capacitance_per_area(tox_nm, k_ox)
    s = _sign_dope(dopant_type)

    # cumulative integrals, both starting at zero on the first grid point
    Q_dep = Q * _cumtrapz0(N, x)                    # C/cm^2
    phi_s = (Q / eps_s) * _cumtrapz0(x * N, x)      # V

    good = x > 0
    x, Q_dep, phi_s = x[good], Q_dep[good], phi_s[good]

    Vg = VFB + s * (phi_s + Q_dep / Cox_pa)
    Cs_pa = eps_s / x
    return Vg, total_capacitance_series_per_area(Cox_pa, Cs_pa) * area_cm2


def _cumtrapz0(y, x):
    """Cumulative trapezoidal integral of y dx, same length as x, first
    element zero. (scipy.integrate.cumulative_trapezoid with initial=0,
    written out so the module keeps working on older scipy.)"""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


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


# ===========================================================================
# PART 8 - the Fermi-Dirac formulation (Seiwatz & Green)
# ---------------------------------------------------------------------------
# Everything above this line uses Boltzmann statistics and assumes the
# dopants are fully ionised. That is the standard textbook treatment and it
# is accurate to a few percent for moderately doped silicon at room
# temperature. This part replaces both assumptions with the degenerate
# treatment: carrier densities from Fermi-Dirac integrals, and dopant
# ionisation from the occupation of the dopant level itself.
#
# The chain is the one set out in the source theory document:
#   1. solve bulk charge neutrality for E_F                  bulk_fermi_level
#   2. assume phi_s, get the surface densities               carrier_densities_fd
#   3. integrate Poisson once for Q_s                        space_charge_density_fd
#   4. C_s = -dQ_s/dphi_s, minority carriers excluded at HF  semiconductor_capacitance_fd
#   5. charge balance gives V_G for that phi_s               gate_voltage_from_phi_s
#   6. sweep phi_s, and the (V_G, C) pairs are the C-V curve cv_curve_fd
# ===========================================================================

# Density-of-states effective masses (m*/m0) and the resulting band-edge
# densities of states at 300 K. N_c = 2*(2*pi*m*kT/h^2)^(3/2).
ME_DOS_SI = 1.09
MH_DOS_SI = 1.15
H_PLANCK = 6.62607015e-34     # J s
M0 = 9.1093837015e-31         # kg

# Shallow dopant ionisation energies in silicon, eV
EA_BORON_EV = 0.045           # above E_v
ED_PHOS_EV = 0.045            # below E_c

# Degeneracy factors for the dopant levels. The source theory document uses
# 2 for both; the more usual textbook choice is 4 for acceptors (the
# valence band is degenerate) and 2 for donors. Both are exposed so the
# difference can be shown rather than argued about - it matters only when
# the dopant level is close to E_F, i.e. at high doping or low temperature.
G_DONOR = 2.0
G_ACCEPTOR = 2.0


def band_dos(T=T_ROOM, m_e=ME_DOS_SI, m_h=MH_DOS_SI):
    """
    Conduction- and valence-band effective densities of states (cm^-3) at
    temperature T, from N = 2*(2*pi*m* k T / h^2)^(3/2).

    Returns (N_c, N_v). At 300 K with the masses above this gives
    N_c ~ 2.8e19 and N_v ~ 3.1e19 cm^-3, the standard silicon values.
    """
    pref = 2.0 * (2.0 * np.pi * KB * T / H_PLANCK**2) ** 1.5
    return pref * (m_e * M0) ** 1.5 * 1e-6, pref * (m_h * M0) ** 1.5 * 1e-6


def fermi_dirac_integral(j, eta, n_quad=96):
    """
    Normalised Fermi-Dirac integral of order j:

        F_j(eta) = (1/Gamma(j+1)) * int_0^inf  e^j / (1 + exp(e - eta)) de

    Normalised so that F_j(eta) -> exp(eta) in the non-degenerate limit,
    which makes n = N_c * F_{1/2}(eta) reduce cleanly to the Boltzmann
    result n = N_c exp((E_F-E_c)/kT).

    The source theory document writes the *unnormalised* integral and
    carries the prefactor 4*pi*(2 m* k T/h^2)^{3/2} explicitly; the two
    forms are related by F_j^unnorm = Gamma(j+1) * F_j^norm, and for
    j = 1/2, 4*pi*(2 m* kT/h^2)^{3/2} * Gamma(3/2) = N_c.

    Evaluated by Gauss-Legendre quadrature on a substitution that maps
    (0, inf) to (0, 1), which is accurate to ~1e-12 over the range of eta
    that matters here (-40 to +40) and vectorises over eta.
    """
    eta = np.asarray(eta, dtype=float)
    scalar = eta.ndim == 0
    eta = np.atleast_1d(eta)
    # substitution e = t/(1-t) maps (0,1) -> (0,inf), de = dt/(1-t)^2
    t, w = np.polynomial.legendre.leggauss(n_quad)
    t = 0.5 * (t + 1.0)
    w = 0.5 * w
    e = t / (1.0 - t)
    jac = 1.0 / (1.0 - t) ** 2
    arg = e[None, :] - eta[..., None]
    with np.errstate(over='ignore'):
        occ = 1.0 / (1.0 + np.exp(np.clip(arg, -700, 700)))
    integrand = e[None, :] ** j * occ * jac[None, :]
    from math import gamma
    out = np.sum(integrand * w[None, :], axis=-1) / gamma(j + 1.0)
    return float(out[0]) if scalar else out


def ionised_dopants(Ef_minus_Ei_eV, NA_cm3, ND_cm3, T=T_ROOM,
                    Eg_eV=EG_SI_EV, EA_eV=EA_BORON_EV, ED_eV=ED_PHOS_EV,
                    g_a=G_ACCEPTOR, g_d=G_DONOR):
    """
    Ionised acceptor and donor densities, allowing for incomplete
    ionisation:

        N_A^- = N_A / (1 + g_a exp((E_A - E_F)/kT))
        N_D^+ = N_D / (1 + g_d exp((E_F - E_D)/kT))

    Energies are referred to the intrinsic level E_i, taken at mid-gap.
    Returns (N_A_ionised, N_D_ionised) in cm^-3.
    """
    Vt = thermal_voltage(T)
    u = np.asarray(Ef_minus_Ei_eV, dtype=float)
    EA = -Eg_eV / 2.0 + EA_eV        # acceptor level, relative to E_i
    ED = Eg_eV / 2.0 - ED_eV         # donor level, relative to E_i
    NA_ion = NA_cm3 / (1.0 + g_a * np.exp(np.clip((EA - u) / Vt, -700, 700)))
    ND_ion = ND_cm3 / (1.0 + g_d * np.exp(np.clip((u - ED) / Vt, -700, 700)))
    return NA_ion, ND_ion


def carrier_densities_fd(Ef_minus_Ei_eV, T=T_ROOM, Eg_eV=EG_SI_EV,
                         Nc=None, Nv=None):
    """
    Electron and hole densities from Fermi-Dirac statistics:

        n = N_c F_{1/2}((E_F - E_c)/kT),  p = N_v F_{1/2}((E_v - E_F)/kT)

    with E_F given relative to the intrinsic level, E_c = E_i + Eg/2 and
    E_v = E_i - Eg/2. Returns (n, p) in cm^-3.
    """
    Vt = thermal_voltage(T)
    if Nc is None or Nv is None:
        Nc, Nv = band_dos(T)
    u = np.asarray(Ef_minus_Ei_eV, dtype=float)
    n = Nc * fermi_dirac_integral(0.5, (u - Eg_eV / 2.0) / Vt)
    p = Nv * fermi_dirac_integral(0.5, (-Eg_eV / 2.0 - u) / Vt)
    return n, p


@lru_cache(maxsize=256)
def _bulk_fermi_level_cached(NA_cm3, ND_cm3, T, Eg_eV, kw_items):
    return _bulk_fermi_level_impl(NA_cm3, ND_cm3, T, Eg_eV, **dict(kw_items))


def bulk_fermi_level(NA_cm3, ND_cm3, T=T_ROOM, Eg_eV=EG_SI_EV, **kw):
    """Cached wrapper - see _bulk_fermi_level_impl for the physics. The
    solve is a brentq over a quadrature and was being repeated at every
    point of every sweep; memoising it cut the notebook runtime by an
    order of magnitude."""
    return _bulk_fermi_level_cached(float(NA_cm3), float(ND_cm3), float(T),
                                    float(Eg_eV), tuple(sorted(kw.items())))


def _bulk_fermi_level_impl(NA_cm3, ND_cm3, T=T_ROOM, Eg_eV=EG_SI_EV, **kw):
    """
    Bulk Fermi level from charge neutrality - step 1 of the procedure:

        N_D^+ - N_A^- + p - n = 0

    solved for E_F by bisection. Returns E_F - E_i in eV.

    This is the step the Boltzmann treatment skips by writing
    phi_F = V_t ln(N/n_i) and assuming full ionisation. At 1e15-1e17 cm^-3
    the two agree closely; the difference grows with doping.
    """
    Nc, Nv = band_dos(T)

    def neutrality(u):
        n, p = carrier_densities_fd(u, T, Eg_eV, Nc, Nv)
        NA_i, ND_i = ionised_dopants(u, NA_cm3, ND_cm3, T, Eg_eV, **kw)
        return ND_i - NA_i + p - n

    return brentq(neutrality, -Eg_eV, Eg_eV, xtol=1e-13, rtol=1e-15)


def intrinsic_density_fd(T=T_ROOM, Eg_eV=EG_SI_EV):
    """
    Intrinsic density implied by the DOS model, n_i = sqrt(N_c N_v)
    exp(-Eg/2kT).

    Worth computing explicitly because it does NOT agree with the
    NI_SI = 9.65e9 cm^-3 used by the Boltzmann half of this module: the
    simple parabolic-band DOS with Eg = 1.12 eV gives about 1.16e10, some
    20% higher. The tabulated value comes from measurement and folds in
    band-structure detail this model does not have. The two halves of the
    module are each internally consistent; do not mix their n_i.
    """
    Nc, Nv = band_dos(T)
    return np.sqrt(Nc * Nv) * np.exp(-Eg_eV / (2 * thermal_voltage(T)))


def _charge_terms_fd(phi, Ef_bulk_eV, NA_cm3, ND_cm3, T, Eg_eV, kw):
    """n, p, N_A^-, N_D^+ at band bending phi (V). Bending phi moves every
    band edge down in energy by q*phi relative to E_F, so it enters as a
    shift of the local (E_F - E_i)."""
    u_local = np.asarray(Ef_bulk_eV, dtype=float) + np.asarray(phi, float)
    n, p = carrier_densities_fd(u_local, T, Eg_eV)
    NA_i, ND_i = ionised_dopants(u_local, NA_cm3, ND_cm3, T, Eg_eV, **kw)
    return n, p, NA_i, ND_i


def space_charge_density_fd(phi_s, NA_cm3, ND_cm3, T=T_ROOM,
                            Eg_eV=EG_SI_EV, eps_s=EPS_SI, n_quad=160,
                            minority_frozen=False, **kw):
    """
    Space-charge density Q_s (C/cm^2) from the Fermi-Dirac formulation -
    step 3 of the procedure.

    Integrating Poisson once from the neutral bulk to the surface,

        xi_s^2 = (2q/eps_s) * int_0^phi_s [ n - p - N_D^+ + N_A^- ] dphi
        Q_s    = -sign(phi_s) * eps_s * |xi_s|

    The integrand is the net charge density that has to be balanced, with
    every term a function of the local band bending. The source theory
    document gives this integral in closed form using F_{3/2} Fermi
    integrals; evaluating it numerically is the same physics and avoids
    transcribing a long closed form, at a cost of one quadrature per point.

    minority_frozen=True drops the minority-carrier term from the
    integrand, which is what the high-frequency measurement sees: minority
    carriers cannot be generated fast enough to follow the AC signal, so
    they contribute no charge to dQ_s/dphi_s.

    Verified against space_charge_density() (the Boltzmann, fully-ionised
    version) in the non-degenerate limit - see cgv_helper_checks.py.
    """
    Ef_b = bulk_fermi_level(NA_cm3, ND_cm3, T, Eg_eV, **kw)
    p_type = NA_cm3 >= ND_cm3
    phi_s = np.atleast_1d(np.asarray(phi_s, dtype=float))

    t, w = _leggauss_cached(n_quad)
    n_b0, p_b0, _, _ = _charge_terms_fd(0.0, Ef_b, NA_cm3, ND_cm3, T, Eg_eV, kw)
    out = np.empty_like(phi_s)
    for i, ps in enumerate(phi_s):
        if ps == 0.0:
            out[i] = 0.0
            continue
        phi = 0.5 * ps * (t + 1.0)                 # 0 .. phi_s
        jac = 0.5 * ps
        n, p, NA_i, ND_i = _charge_terms_fd(phi, Ef_b, NA_cm3, ND_cm3,
                                            T, Eg_eV, kw)
        if minority_frozen:
            if p_type:
                n = np.full_like(n, float(np.asarray(n_b0).ravel()[0]))
            else:
                p = np.full_like(p, float(np.asarray(p_b0).ravel()[0]))
        integ = np.sum((n - p - ND_i + NA_i) * w) * jac
        xi2 = 2.0 * Q / eps_s * integ
        out[i] = -np.sign(ps) * eps_s * np.sqrt(max(xi2, 0.0))
    return out


def semiconductor_capacitance_fd(phi_s, NA_cm3, ND_cm3, T=T_ROOM,
                                 Eg_eV=EG_SI_EV, h=2e-4,
                                 minority_frozen=False, **kw):
    """
    C_s = -dQ_s/dphi_s (F/cm^2) from the Fermi-Dirac Q_s, by central
    difference. Same minus sign as semiconductor_capacitance(): Q_s falls
    as phi_s rises, so the bare derivative is negative.
    """
    phi_s = np.atleast_1d(np.asarray(phi_s, dtype=float))
    kwargs = dict(T=T, Eg_eV=Eg_eV, minority_frozen=minority_frozen, **kw)
    qp = space_charge_density_fd(phi_s + h, NA_cm3, ND_cm3, **kwargs)
    qm = space_charge_density_fd(phi_s - h, NA_cm3, ND_cm3, **kwargs)
    return -(qp - qm) / (2 * h)


def interface_charge_fd(phi_s, Dit_donor, Dit_acceptor, NA_cm3, ND_cm3,
                        sigma_ratio=1.0, T=T_ROOM, Eg_eV=EG_SI_EV,
                        n_points=201, **kw):
    """
    Interface-trap charge Q_it (C/cm^2) with Shockley-Read-Hall occupancy:

        f_d = [(sn/sp) n1 + p_s] / [(sn/sp)(n_s + n1) + (p_s + p1)]
        Q_it = q * int D_it,donor f_d dE  -  q * int D_it,acceptor (1-f_d) dE

    Donor-like states are positive when empty of electrons, acceptor-like
    negative when full, and both occupancies are set by the *surface*
    carrier densities, so Q_it moves as the bands bend. This is the same
    statistics used in KPSPV/kpspv_helper.py, and deliberately so - the two
    techniques are looking at the same interface.

    Dit_donor and Dit_acceptor are callables D(E) taking energy above E_v
    in eV and returning cm^-2 eV^-1.
    """
    Vt = thermal_voltage(T)
    Ef_b = bulk_fermi_level(NA_cm3, ND_cm3, T, Eg_eV, **kw)
    ni = intrinsic_density_fd(T, Eg_eV)
    phi_s = np.atleast_1d(np.asarray(phi_s, dtype=float))

    E = np.linspace(1e-4, Eg_eV - 1e-4, n_points)      # from E_v
    Ei = Eg_eV / 2.0
    n1 = ni * np.exp(np.clip((E - Ei) / Vt, -700, 700))
    p1 = ni * np.exp(np.clip((Ei - E) / Vt, -700, 700))
    Dd = np.asarray(Dit_donor(E), dtype=float)
    Da = np.asarray(Dit_acceptor(E), dtype=float)

    out = np.empty_like(phi_s)
    for i, ps in enumerate(phi_s):
        n_s, p_s, _, _ = _charge_terms_fd(ps, Ef_b, NA_cm3, ND_cm3,
                                          T, Eg_eV, kw)
        f_d = (sigma_ratio * n1 + p_s) / (sigma_ratio * (n_s + n1)
                                          + (p_s + p1))
        out[i] = Q * (_trapz(Dd * f_d, E) - _trapz(Da * (1.0 - f_d), E))
    return out


def gate_voltage_from_phi_s(phi_s, Cox_pa, phi_ms_eV, Qs_C_cm2,
                            Qit_C_cm2=0.0, Qf_C_cm2=0.0, centroid_cm=0.0,
                            tox_cm=None):
    """
    Gate voltage for an assumed surface potential - step 5 of the
    procedure, and the equation the whole C-V calculation hangs on:

        V_G = phi_ms/q + phi_s - [ Q_s + Q_it + (1 + d/t_i) Q_f ] / C_i

    The (1 + d/t_i) factor is the charge-centroid weighting. d is measured
    from the insulator-semiconductor interface towards the gate and runs
    from 0 to -t_i, so the weight is 1 for charge sitting right at the
    semiconductor interface (full effect on V_G) and 0 for charge at the
    gate (no effect at all). Fixed charge only shifts the curve to the
    extent that it is separated from the semiconductor.

    Pass centroid_cm and tox_cm to use the weighting; leaving them at the
    default treats Q_f as sitting at the interface (weight 1), which is
    the usual textbook simplification.
    """
    weight = 1.0
    if tox_cm:
        weight = 1.0 + np.asarray(centroid_cm, dtype=float) / tox_cm
    return (phi_ms_eV + np.asarray(phi_s, dtype=float)
            - (np.asarray(Qs_C_cm2, dtype=float)
               + np.asarray(Qit_C_cm2, dtype=float)
               + weight * np.asarray(Qf_C_cm2, dtype=float)) / Cox_pa)


def cv_curve_fd(phi_s_grid, area_cm2, tox_nm, NA_cm3, ND_cm3,
                phi_ms_eV=0.0, Qf_C_cm2=0.0, centroid_cm=0.0,
                Dit_donor=None, Dit_acceptor=None, sigma_ratio=1.0,
                high_frequency=True, T=T_ROOM, Eg_eV=EG_SI_EV,
                k_ox=K_SIO2, **kw):
    """
    The complete theoretical C-V curve, following the procedure in the
    source theory document end to end - steps 1 to 6.

    Rather than solving for phi_s at each V_G (what surface_potential_from_bias
    does for the Boltzmann model), this sweeps phi_s and computes the V_G
    that produces it. That is the natural direction: every quantity is an
    explicit function of phi_s, so no root-finding is needed and the curve
    comes out parametrically.

    Returns (V_G, C_total, phi_s, Q_s, Q_it), all arrays.

    high_frequency=True freezes the minority carriers in dQ_s/dphi_s (the
    HF curve). The DC charge balance that sets V_G always uses the full
    Q_s including minority carriers, because that is a DC quantity - this
    asymmetry is the whole content of the HF/LF distinction and is easy to
    get wrong by freezing them in both places.
    """
    Cox_pa = oxide_capacitance_per_area(tox_nm, k_ox)
    tox_cm = tox_nm * 1e-7
    phi = np.atleast_1d(np.asarray(phi_s_grid, dtype=float))

    Qs_dc = space_charge_density_fd(phi, NA_cm3, ND_cm3, T, Eg_eV,
                                    minority_frozen=False, **kw)
    Cs = semiconductor_capacitance_fd(phi, NA_cm3, ND_cm3, T, Eg_eV,
                                      minority_frozen=high_frequency, **kw)
    if Dit_donor is not None and Dit_acceptor is not None:
        Qit = interface_charge_fd(phi, Dit_donor, Dit_acceptor, NA_cm3,
                                  ND_cm3, sigma_ratio, T, Eg_eV, **kw)
    else:
        Qit = np.zeros_like(phi)

    Vg = gate_voltage_from_phi_s(phi, Cox_pa, phi_ms_eV, Qs_dc, Qit,
                                 Qf_C_cm2, centroid_cm, tox_cm)
    C = total_capacitance_series_per_area(Cox_pa, Cs) * area_cm2
    return Vg, C, phi, Qs_dc, Qit


@lru_cache(maxsize=32)
def _phi_qs_table(NA_cm3, ND_cm3, T, Eg_eV, kw_items, n=400):
    """Cached phi_s grid and the matching Q_s, used to invert the charge
    balance by interpolation rather than by a root find per sample."""
    phi = np.linspace(-Eg_eV - 0.4, Eg_eV + 0.4, n)
    return phi, space_charge_density_fd(phi, NA_cm3, ND_cm3, T, Eg_eV,
                                        **dict(kw_items))


def phi_s_distribution(Vg, Cox_pa, phi_ms_eV, Qf0_C_cm2, sigma_q_C_cm2,
                       NA_cm3, ND_cm3, n_sigma=3.0, n_points=21,
                       T=T_ROOM, Eg_eV=EG_SI_EV, **kw):
    """
    Distribution of surface potential produced by a fluctuating fixed
    charge - the "how to integrate charge fluctuations" step of the source
    theory document.

    Real dielectric films are not uniformly charged. Writing
    Q_f' = Q_f + dq and re-solving

        V_G - phi_ms/q - Q_f'/C_i = phi_s' - Q_s(phi_s')/C_i

    for each sample of a Gaussian distribution in Q_f gives a distribution
    P(phi_s) rather than a single value. A large-area gate then measures
    the average over that distribution, which smears every feature of the
    C-V curve exactly as it smears the KP charge sweep in the KP/SPV
    notebook (Eq. 21 there).

    Returns (phi_s_samples, weights) with the weights normalised to 1.
    """
    q_grid = np.linspace(Qf0_C_cm2 - n_sigma * sigma_q_C_cm2,
                         Qf0_C_cm2 + n_sigma * sigma_q_C_cm2, n_points)
    w = np.exp(-0.5 * ((q_grid - Qf0_C_cm2) / sigma_q_C_cm2) ** 2)
    w = w / w.sum()

    # Build the monotonic map phi_s -> V_G once and invert it by
    # interpolation, instead of running a root find (each of which calls
    # the Q_s quadrature dozens of times) for every charge sample. The map
    # is strictly increasing in phi_s, so np.interp is exact to the grid.
    phi_tab, Qs_tab = _phi_qs_table(float(NA_cm3), float(ND_cm3), float(T),
                                    float(Eg_eV), tuple(sorted(kw.items())))
    base = phi_ms_eV + phi_tab - Qs_tab / Cox_pa      # V_G at Q_f = 0
    order = np.argsort(base)
    out = np.interp(Vg + q_grid / Cox_pa, base[order], phi_tab[order])
    return out, w


# ---------------------------------------------------------------------------
# The full interface-state admittance: every branch of the equivalent
# circuit, integrated over the gap. Sec. 13's Lorentzian is the
# single-level limit of this.
# ---------------------------------------------------------------------------
def interface_branch_admittance(omega, phi_s, Dit_func, NA_cm3, ND_cm3,
                                sigma_n_cm2=1e-16, sigma_p_cm2=1e-16,
                                v_th=VTH_SI, T=T_ROOM, Eg_eV=EG_SI_EV,
                                n_points=201, **kw):
    """
    Capacitance and conductance of the interface-state network, summed
    over the band gap, following the source theory document's Sec. 2.

    For a trap at energy E with occupancy f_t, the two capture paths -
    to the valence band and to the conduction band - have resistances

        R_ps = V_t / (q f_t S_p0 p_s),   S_p0 = sigma_p v_th D_it
        R_ns = V_t / (q (1-f_t) S_n0 n_s), S_n0 = sigma_n v_th D_it

    and the trap itself a capacitance C_it = q D_it f_t (1-f_t)/V_t.
    Those three elements give six frequency-dependent quantities
    (C_dp, C_dn, C_pn and G_dp, G_dn, G_pn) once integrated over E.

    Returns a dict with those six arrays, shaped like omega.

    The single-level Lorentzian of Eq. (17) is what this collapses to when
    D_it is a delta function and one capture path dominates. The point of
    carrying all of it is that a real interface has a *distribution* of
    time constants, which is why measured conductance peaks are broader
    and lower than the single-level formula predicts (Sec. 13's second
    curve shows that empirically; this computes it).
    """
    Vt = thermal_voltage(T)
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    Ef_b = bulk_fermi_level(NA_cm3, ND_cm3, T, Eg_eV, **kw)
    n_s, p_s, _, _ = _charge_terms_fd(phi_s, Ef_b, NA_cm3, ND_cm3, T,
                                      Eg_eV, kw)
    n_s, p_s = float(np.asarray(n_s).ravel()[0]), float(np.asarray(p_s).ravel()[0])

    E = np.linspace(1e-3, Eg_eV - 1e-3, n_points)     # from E_v
    Dit = np.asarray(Dit_func(E), dtype=float)
    # trap occupancy referred to the surface quasi-Fermi level
    Ef_surface = Ef_b + float(np.asarray(phi_s).ravel()[0]) + Eg_eV / 2.0
    f_t = 1.0 / (1.0 + np.exp(np.clip((E - Ef_surface) / Vt, -700, 700)))

    Sp0 = sigma_p_cm2 * v_th * Dit
    Sn0 = sigma_n_cm2 * v_th * Dit
    with np.errstate(divide='ignore', invalid='ignore'):
        R_ps = Vt / (Q * np.clip(f_t, 1e-300, None) * Sp0 * p_s)
        R_ns = Vt / (Q * np.clip(1 - f_t, 1e-300, None) * Sn0 * n_s)
    C_it = Q * Dit * f_t * (1.0 - f_t) / Vt

    w = omega[:, None]
    r = R_ps / R_ns
    den_p = (1.0 + r) ** 2 + (w * R_ps * C_it) ** 2
    den_n = (1.0 + 1.0 / r) ** 2 + (w * R_ns * C_it) ** 2
    den_pn = (R_ps + R_ns) ** 2 + (w * R_ps * R_ns * C_it) ** 2

    out = {
        'C_dp': _trapz((1.0 + r) * C_it / den_p, E, axis=-1),
        'C_dn': _trapz((1.0 + 1.0 / r) * C_it / den_n, E, axis=-1),
        'C_pn': _trapz(R_ns * R_ps * C_it / den_pn, E, axis=-1),
        'G_dp': _trapz(w**2 * R_ps * C_it**2 / den_p, E, axis=-1),
        'G_dn': _trapz(w**2 * R_ns * C_it**2 / den_n, E, axis=-1),
        'G_pn': _trapz((R_ps + R_ns) / den_pn, E, axis=-1),
    }
    return out


def interface_conductance_parallel(branches):
    """
    The trap-response conductance, G_p = G_dp + G_dn (S/cm^2).

    These are the two branches in which a trap exchanges carriers with a
    band and dissipates energy doing it - the loss the conductance method
    is built to measure. G_p/omega peaks at omega*tau = 1 for the branch
    that dominates, which on a depleted p-type surface is the
    trap-to-valence-band path with tau_p = 1/(sigma_p v_th p_s). Because
    p_s is set by the band bending, the peak scans in frequency as the
    bias is swept, and that is the whole technique.

    Use this, not total_interface_conductance(), to reproduce a measured
    G_p/omega peak - see the note on that function.
    """
    return branches['G_dp'] + branches['G_dn']


def total_interface_conductance(omega, Cox_pa, branches):
    """
    Total conductance of the full six-element network, transcribed from the
    last equation of the source theory document's conductance section:

        G_mT = G_dn + [ (G_dp+G_pn)(G_dp G_pn - w^2 C_pn (C_i+C_dp))
                        + w^2 (C_i+C_dp+C_pn)(G_dp C_pn + G_pn (C_i+C_dp)) ]
                      / [ (G_dp+G_pn)^2 + w^2 (C_i+C_dp+C_pn)^2 ]

    `branches` is the dict returned by interface_branch_admittance.

    CAUTION - this does not reproduce a conductance peak, and it is worth
    knowing why before using it. The expression includes G_pn, the
    band-to-band path *through* the traps, which is a DC recombination
    conductance: it does not vanish as omega -> 0. So G/omega diverges at
    low frequency and completely masks the trap-response peak that the
    conductance method is looking for. Numerically, at phi_s = 0.2 V on
    1e16 cm^-3 p-type with D_it = 1e11, G_dp/omega peaks cleanly at
    1.5 kHz, while this total falls monotonically from 1 Hz upwards.

    That is a real feature of the network, not a coding error - but it
    means this expression answers a different question from the one the
    conductance method asks. For the measured peak use
    interface_conductance_parallel(); this is here for completeness and
    for the case where band-to-band recombination through the interface is
    itself the quantity of interest.
    """
    w = np.asarray(omega, dtype=float)
    Gdp, Gdn, Gpn = branches['G_dp'], branches['G_dn'], branches['G_pn']
    Cdp, Cpn = branches['C_dp'], branches['C_pn']
    A = Cox_pa + Cdp
    num = ((Gdp + Gpn) * (Gdp * Gpn - w**2 * Cpn * A)
           + w**2 * (A + Cpn) * (Gdp * Cpn + Gpn * A))
    den = (Gdp + Gpn) ** 2 + w**2 * (A + Cpn) ** 2
    return Gdn + num / den
