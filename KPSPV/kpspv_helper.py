"""
kpspv_helper.py
---------------
Helper functions for the Kelvin probe / surface photovoltage (KP/SPV)
teaching notebook, `kpspv_analysis.ipynb`. Equation numbers quoted in the
docstrings below, e.g. "Eq. (12)", refer to the numbered equations in that
notebook.

Scope
-----
These functions implement a 1-D electrostatic description of the system a
macro-scale Kelvin probe actually looks at: a vibrating metal electrode
above a semiconductor that may carry a dielectric film, with charge stored
in the film, in interface states, and in the semiconductor space-charge
region. Illumination enters through a single excess-carrier density.

They are written for teaching, not for precision surface metrology:
  - the surface is uniform and one-dimensional, while a millimetre-scale
    probe averages over a patchy surface,
  - carrier statistics are Boltzmann, not Fermi-Dirac, and the dopants are
    fully ionised,
  - the excess carrier density is taken as constant across the space-charge
    region,
  - band gap narrowing, surface dipoles and image-force lowering are all
    ignored,
  - the interface-state distribution is a smooth model shape, not a
    measured spectrum.

For quantitative work, use tabulated silicon parameters, a Fermi-Dirac
solver, a measured D_it spectrum, and an optical model for the generation
profile.

SIGN CONVENTION - read this before using anything below
-------------------------------------------------------
Both CPD and SPV appear with either sign in the literature. This module
uses one convention throughout:

  * The Kelvin probe applies a backing voltage V_b to the probe and finds
    the value that nulls the signal. The contact potential difference is
        CPD = -V_b(null)
    so for two bare conductors  CPD = (Phi_sample - Phi_probe) / q.

  * The semiconductor surface potential phi_s is measured relative to the
    neutral bulk, with the surface carrier densities
        n_s = n_b exp(+phi_s / V_t),   p_s = p_b exp(-phi_s / V_t).
    On n-type material phi_s > 0 therefore means ACCUMULATION and
    phi_s < 0 means depletion or inversion. (Some reviews define phi_s with
    the opposite sign on p-type material - do not mix the two.)

  * The surface photovoltage is
        SPV = CPD_dark - CPD_light
    which, when the light drives the surface to flat band, equals
    -phi_s(dark).

UNIT CONVENTION
---------------
Energies and work functions are in eV. Potentials and voltages are in V.
Carrier and dopant densities are in cm^-3, sheet charge densities in
elementary charges per cm^2 (signed, so -2e11 means 2e11 electrons per
cm^2). Dielectric thicknesses and charge centroids are in nm, wafer
thickness and depletion widths in um, probe spacings in um, probe area in
mm^2. Conversions happen inside these functions, never in the notebook.
"""

import numpy as np
from scipy.optimize import brentq

# numpy >= 2.0 renamed trapz -> trapezoid; keep this notebook working on
# either version.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ---------------------------------------------------------------------------
# Physical constants (cm-based electrostatics, as is conventional for
# semiconductor surface work)
# ---------------------------------------------------------------------------
Q = 1.602176634e-19       # elementary charge, C
KB = 1.380649e-23         # Boltzmann constant, J/K
EPS0 = 8.8541878128e-14   # vacuum permittivity, F/cm

T_ROOM = 300.0            # K

# Crystalline silicon at 300 K
EG_SI_EV = 1.12           # band gap, eV
CHI_SI_EV = 4.05          # electron affinity, eV (reported 4.05-4.15;
                          # this spread alone is +/-50 meV on any absolute
                          # semiconductor work function computed from Eq. 8)
K_SI = 11.7               # relative permittivity
NI_SI = 9.65e9            # intrinsic carrier density, cm^-3
VTH_SI = 1.1e7            # thermal velocity, cm/s

# Common dielectrics
K_SIO2 = 3.9
K_SIN = 7.5
K_AL2O3 = 9.0

# Reference work functions, eV. The spread between reported values is of
# order 100 meV - these are nominal figures for calibration examples only.
WORK_FUNCTION_EV = {
    'Au': 5.10,
    'Pt': 5.65,
    'Ni': 5.15,
    'Cu': 4.65,
    'W': 4.55,
    'stainless steel': 4.40,
    'Al': 4.28,
    'HOPG': 4.60,
}


def thermal_voltage(T=T_ROOM):
    """Thermal voltage V_t = kT/q, in volts."""
    return KB * T / Q


# ===========================================================================
# PART 1 - the Kelvin probe itself
# ===========================================================================
def cpd_from_work_functions(phi_sample_eV, phi_probe_eV):
    """
    Contact potential difference between two bare conductors - Eq. (2).

    CPD = (Phi_sample - Phi_probe) / q, in volts, with both work functions
    in eV. Positive CPD means the sample has the higher work function.
    """
    return np.asarray(phi_sample_eV, dtype=float) - float(phi_probe_eV)


def work_function_from_cpd(cpd_V, phi_probe_eV):
    """
    Absolute sample work function from a measured CPD - Eq. (7).

    Phi_sample = Phi_probe + q*CPD. This is the whole reason a Kelvin probe
    needs calibrating: the measurement gives a difference, and the probe
    work function has to come from somewhere else.
    """
    return float(phi_probe_eV) + np.asarray(cpd_V, dtype=float)


def calibrate_probe(cpd_reference_V, phi_reference_eV):
    """
    Probe work function from a measurement on a reference sample - Eq. (7)
    rearranged: Phi_probe = Phi_reference - q*CPD_reference.
    """
    return float(phi_reference_eV) - float(cpd_reference_V)


def modulation_index(d0_um, d1_um):
    """
    Modulation index eps = d1/d0 of the vibrating capacitor: the peak
    excursion as a fraction of the mean spacing. Signal grows with eps, and
    so does the harmonic content of the current waveform.
    """
    return float(d1_um) / float(d0_um)


def kelvin_capacitance(t_s, area_mm2=3.14, d0_um=200.0, d1_um=40.0,
                       freq_Hz=80.0):
    """
    Capacitance of the vibrating Kelvin capacitor - Eq. (3).

        C(t) = eps0 * A / (d0 + d1 sin(wt))

    Returns capacitance in farads for time(s) in seconds. Defaults are a
    2 mm diameter probe tip at a 200 um mean spacing - a typical ambient
    macro-scale probe.
    """
    t_s = np.asarray(t_s, dtype=float)
    area_cm2 = area_mm2 * 1e-2                 # mm^2 -> cm^2
    d_cm = (d0_um + d1_um * np.sin(2 * np.pi * freq_Hz * t_s)) * 1e-4
    return EPS0 * area_cm2 / d_cm


def kelvin_current(t_s, V_backing, cpd_V, area_mm2=3.14, d0_um=200.0,
                   d1_um=40.0, freq_Hz=80.0):
    """
    Current through the Kelvin circuit - Eq. (4).

        i(t) = (V_b + CPD) * dC/dt

    The DC quantity CPD is converted into an AC current purely by the
    mechanical modulation of C. The current vanishes for all t when
    V_b = -CPD, which is the null condition of Eq. (5).
    """
    t_s = np.asarray(t_s, dtype=float)
    w = 2 * np.pi * freq_Hz
    area_cm2 = area_mm2 * 1e-2
    d_cm = (d0_um + d1_um * np.sin(w * t_s)) * 1e-4
    dd_dt = d1_um * 1e-4 * w * np.cos(w * t_s)
    dC_dt = -EPS0 * area_cm2 * dd_dt / d_cm ** 2
    return (np.asarray(V_backing, dtype=float) + cpd_V) * dC_dt


def offnull_amplitude(V_backing, cpd_V, gradient_V_per_V=3.0):
    """
    Peak-to-peak Kelvin signal against backing voltage - Eq. (5).

        V_ptp = k * (V_b + CPD)

    Linear in V_b, with a zero crossing at V_b = -CPD and a phase inversion
    through it. The gradient k depends on the probe area, the mean spacing
    and the amplifier gain - all of which drift - which is exactly why the
    zero crossing, and not the amplitude, is the measured quantity.
    """
    return gradient_V_per_V * (np.asarray(V_backing, dtype=float) + cpd_V)


class NullFitResult(dict):
    """
    Result of an off-null Kelvin fit. A plain dict with attribute access:
    `res.cpd_V` and `res['cpd_V']` both work.

    Keys: cpd_V, cpd_err, gradient, gradient_err, r_squared, V_backing,
    V_ptp.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def summary(self):
        """A formatted report of the extracted contact potential."""
        return "\n".join([
            "Off-null Kelvin fit",
            "-------------------",
            f"  fit        V_ptp = {self['gradient']:.4g} * V_b "
            f"+ {self['gradient'] * self['cpd_V']:.4g}"
            f"   (R^2 = {self['r_squared']:.5f})",
            f"  gradient   k     = {self['gradient']:.4g} "
            f"+/- {self['gradient_err']:.2g} V/V",
            f"  contact p. CPD   = {self['cpd_V'] * 1e3:.5g} "
            f"+/- {self['cpd_err'] * 1e3:.2g} mV",
        ])

    def __str__(self):
        return self.summary()


def fit_null_point(V_backing, V_ptp, sigma=None):
    """
    The off-null extraction - Eq. (6).

    Fits V_ptp = k*(V_b + CPD) = k*V_b + c to points measured *away* from
    the null, where the signal-to-noise ratio is good, and reports the zero
    crossing CPD = c/k.

    The uncertainty on CPD uses the full fit covariance, not the separate
    standard errors: slope and intercept of a line fitted over a limited
    range are strongly correlated, and CPD is their ratio.

    Returns a NullFitResult.
    """
    x = np.asarray(V_backing, dtype=float)
    y = np.asarray(V_ptp, dtype=float)
    w = None if sigma is None else 1.0 / np.asarray(sigma, dtype=float)

    (k, c), cov = np.polyfit(x, y, 1, w=w, cov=True)
    cpd_V = c / k                                # since c = k * CPD

    # var(CPD) for CPD = c/k, first-order, including cov(k, c)
    d_dk = -c / k ** 2
    d_dc = 1.0 / k
    var = (d_dk ** 2 * cov[0, 0] + d_dc ** 2 * cov[1, 1]
           + 2 * d_dk * d_dc * cov[0, 1])

    resid = y - (k * x + c)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan

    return NullFitResult(cpd_V=cpd_V, cpd_err=np.sqrt(max(var, 0.0)),
                         gradient=k, gradient_err=np.sqrt(cov[0, 0]),
                         r_squared=r2, V_backing=x, V_ptp=y)


# ===========================================================================
# PART 2 - the semiconductor surface
# ===========================================================================
def fermi_potential(doping_cm3, dopant_type='n', T=T_ROOM, ni=NI_SI):
    """
    Bulk Fermi potential relative to the intrinsic level - Eq. (9).

        phi_F = V_t * ln(N_dop / n_i)

    Returned positive for n-type and negative for p-type, following the
    convention of Eq. (8).
    """
    phi = thermal_voltage(T) * np.log(np.asarray(doping_cm3, float) / ni)
    return phi if dopant_type.lower().startswith('n') else -phi


def bulk_carrier_densities(doping_cm3, dopant_type='n', delta_n_cm3=0.0,
                           ni=NI_SI):
    """
    Equilibrium-plus-injection carrier densities in the neutral bulk.

    Returns (n_b, p_b) in cm^-3. Complete ionisation is assumed, and the
    excess density delta_n is added equally to both, so the net ionised
    charge n_b - p_b is unchanged by illumination.
    """
    N = np.asarray(doping_cm3, dtype=float)
    dn = np.asarray(delta_n_cm3, dtype=float)
    if dopant_type.lower().startswith('n'):
        n0 = 0.5 * (N + np.sqrt(N ** 2 + 4 * ni ** 2))
        p0 = ni ** 2 / n0
    else:
        p0 = 0.5 * (N + np.sqrt(N ** 2 + 4 * ni ** 2))
        n0 = ni ** 2 / p0
    return n0 + dn, p0 + dn


def work_function_semiconductor(doping_cm3, dopant_type='n',
                                chi_eV=CHI_SI_EV, Eg_eV=EG_SI_EV, T=T_ROOM,
                                ni=NI_SI):
    """
    Flat-band work function of a uniformly doped semiconductor - Eq. (8).

        Phi_s = chi + Eg/2 - q*phi_F

    chi is the electron affinity at flat band, so the doping term and the
    band-bending term of Eq. (10) are added separately and are not already
    contained in chi.
    """
    return chi_eV + Eg_eV / 2.0 - fermi_potential(doping_cm3, dopant_type,
                                                  T, ni)


def surface_carrier_densities(phi_s_V, n_b, p_b, T=T_ROOM):
    """
    Carrier densities at the semiconductor surface - Eq. (13).

        n_s = n_b exp(+phi_s/V_t),   p_s = p_b exp(-phi_s/V_t)

    On n-type material phi_s > 0 raises n_s (accumulation) and phi_s < 0
    raises p_s (depletion, then inversion).
    """
    u = np.asarray(phi_s_V, dtype=float) / thermal_voltage(T)
    u = np.clip(u, -200.0, 200.0)          # keep exp() finite
    return n_b * np.exp(u), p_b * np.exp(-u)


def space_charge_density(phi_s_V, doping_cm3, dopant_type='n',
                         delta_n_cm3=0.0, K_s=K_SI, T=T_ROOM, ni=NI_SI):
    """
    Charge per unit area in the semiconductor space-charge region - Eq. (11).

        Q_sc = -sign(phi_s) * sqrt(2 q eps_s V_t * G(phi_s)) / q
        G    = n_b(e^u - 1 - u) + p_b(e^-u - 1 + u),   u = phi_s/V_t

    Returned in elementary charges per cm^2, signed: negative means the
    semiconductor holds net negative charge. G is positive definite, so the
    sign of Q_sc is always opposite to the sign of phi_s.

    Accumulation, depletion and inversion, in the dark and under
    illumination, all come out of this one expression - delta_n enters
    through n_b and p_b.
    """
    Vt = thermal_voltage(T)
    n_b, p_b = bulk_carrier_densities(doping_cm3, dopant_type, delta_n_cm3,
                                      ni)
    u = np.clip(np.asarray(phi_s_V, dtype=float) / Vt, -200.0, 200.0)
    G = n_b * (np.expm1(u) - u) + p_b * (np.expm1(-u) + u)
    G = np.maximum(G, 0.0)
    Q_C = -np.sign(u) * np.sqrt(2 * Q * K_s * EPS0 * Vt * G)
    return Q_C / Q                                # C/cm^2 -> charges/cm^2


def band_bending_depletion(sheet_charge_cm2, doping_cm3, K_s=K_SI):
    """
    Surface potential in the depletion approximation - Eq. (12).

        |phi_s| = q N^2 / (2 K_s eps0 N_dop)

    with N the surface charge density in cm^-2. Numerically this is
    9.05e-7 * N^2 / (K_s * N_dop) volts. Useful as a closed-form check on
    Eq. (11), and as the estimate a student can do by hand - but it is
    valid only in depletion, and it fails in accumulation and in strong
    inversion where free carriers, not ionised dopants, hold the charge.
    """
    N = np.abs(np.asarray(sheet_charge_cm2, dtype=float))
    return Q * N ** 2 / (2 * K_s * EPS0 * np.asarray(doping_cm3, float))


def depletion_width_um(phi_s_V, doping_cm3, K_s=K_SI):
    """
    Depletion width for a given band bending, in micrometres:
    W = sqrt(2 K_s eps0 |phi_s| / (q N_dop)). Depletion approximation only.
    """
    phi = np.abs(np.asarray(phi_s_V, dtype=float))
    W_cm = np.sqrt(2 * K_s * EPS0 * phi / (Q * np.asarray(doping_cm3, float)))
    return W_cm * 1e4


def insulator_potential(Qf_cm2, xc_nm, K_i=K_SIO2):
    """
    Potential dropped across a dielectric film by charge inside it - Eq. (14).

        V_i = x_c * Q_f / (K_i eps0)

    Q_f is the effective sheet charge in elementary charges per cm^2 and
    x_c its centroid, measured in nm from the dielectric/semiconductor
    interface. Only the product Q_f * x_c is observable: a Kelvin probe
    cannot separate a small charge far from the interface from a large one
    close to it.
    """
    xc_cm = np.asarray(xc_nm, dtype=float) * 1e-7
    return Q * np.asarray(Qf_cm2, dtype=float) * xc_cm / (K_i * EPS0)


# --- interface states ------------------------------------------------------
def dit_profile(E_eV, dit_midgap=1e11, dit_edge=1e13, Eg_eV=EG_SI_EV,
                tail_eV=0.15):
    """
    Model interface-state density D_it(E), in cm^-2 eV^-1, for energies E
    measured from the valence band edge.

    A flat mid-gap density with exponential band tails rising towards both
    band edges - the U shape reported for the Si/SiO2 interface. A model
    shape, not a measured spectrum.
    """
    E = np.asarray(E_eV, dtype=float)
    if dit_midgap <= 0:
        return np.zeros_like(E)
    ratio = max(dit_edge / dit_midgap, 1.0)
    tails = ratio ** np.exp(-np.minimum(E, Eg_eV - E) / tail_eV)
    return dit_midgap * tails


def _srh_energy_grid(Eg_eV, n_points):
    return np.linspace(1e-4, Eg_eV - 1e-4, n_points)


def interface_charge(phi_s_V, doping_cm3, dopant_type='n', delta_n_cm3=0.0,
                     dit_midgap=1e11, dit_edge=1e13, sigma_n=1e-16,
                     sigma_p=1e-16, Eg_eV=EG_SI_EV, T=T_ROOM, ni=NI_SI,
                     n_points=201):
    """
    Charge stored in interface states - Eqs. (15) and (16).

        Q_it = int D_it,donor * f_p dE  -  int D_it,acceptor * f_n dE

    Donor-like states occupy the lower half of the gap and are positive
    when empty of electrons; acceptor-like states occupy the upper half and
    are negative when full. Occupancy follows Shockley-Read-Hall statistics
    evaluated with the *surface* carrier densities, so Q_it moves as the
    bands bend. That feedback is what makes the CPD-against-charge curve of
    Eq. (17) stretch out along the charge axis.

    Returned in elementary charges per cm^2, signed.
    """
    Vt = thermal_voltage(T)
    n_b, p_b = bulk_carrier_densities(doping_cm3, dopant_type, delta_n_cm3,
                                      ni)
    n_s, p_s = surface_carrier_densities(phi_s_V, n_b, p_b, T)

    E = _srh_energy_grid(Eg_eV, n_points)             # from E_v
    Ei = Eg_eV / 2.0
    n1 = ni * np.exp(np.clip((E - Ei) / Vt, -200, 200))
    p1 = ni * np.exp(np.clip((Ei - E) / Vt, -200, 200))

    r = sigma_n / sigma_p
    f_p = (r * n1 + p_s) / (r * (n_s + n1) + (p_s + p1))   # empty of e-
    f_n = 1.0 - f_p

    dit = dit_profile(E, dit_midgap, dit_edge, Eg_eV)
    donor = np.where(E < Ei, dit, 0.0)
    acceptor = np.where(E >= Ei, dit, 0.0)

    return _trapz(donor * f_p, E) - _trapz(acceptor * f_n, E)


def surface_recombination_velocity(phi_s_V, doping_cm3, dopant_type='n',
                                   delta_n_cm3=1e14, dit_midgap=1e11,
                                   dit_edge=1e13, sigma_n=1e-16,
                                   sigma_p=1e-16, Eg_eV=EG_SI_EV, T=T_ROOM,
                                   ni=NI_SI, n_points=201):
    """
    Surface recombination velocity from the interface-state distribution -
    Eq. (21).

        S = (1/dn) * int (n_s p_s - n_i^2) /
                     [ (n_s + n1)/S_p0 + (p_s + p1)/S_n0 ] dE

    with S_n0,p0 = v_th * sigma_n,p * D_it(E). S in cm/s.

    This is why band bending and passivation are the same subject: pushing
    phi_s away from zero drives one of the two surface carrier densities
    down, the recombination rate collapses with it, and S falls. That is
    field-effect passivation.
    """
    Vt = thermal_voltage(T)
    n_b, p_b = bulk_carrier_densities(doping_cm3, dopant_type, delta_n_cm3,
                                      ni)
    n_s, p_s = surface_carrier_densities(phi_s_V, n_b, p_b, T)

    E = _srh_energy_grid(Eg_eV, n_points)
    Ei = Eg_eV / 2.0
    n1 = ni * np.exp(np.clip((E - Ei) / Vt, -200, 200))
    p1 = ni * np.exp(np.clip((Ei - E) / Vt, -200, 200))

    dit = dit_profile(E, dit_midgap, dit_edge, Eg_eV)
    S_n0 = VTH_SI * sigma_n * dit
    S_p0 = VTH_SI * sigma_p * dit

    integrand = (n_s * p_s - ni ** 2) / ((n_s + n1) / S_p0
                                         + (p_s + p1) / S_n0)
    return _trapz(integrand, E) / delta_n_cm3


def effective_lifetime_us(tau_bulk_us, S_cm_s, W_um):
    """
    Effective lifetime of a wafer passivated on both faces - Eq. (22).

        1/tau_eff = 1/tau_bulk + 2S/W

    Reduces to tau_bulk as S -> 0 and is dominated by the surface term as
    S -> infinity. The same decomposition as the effective diffusion length
    in the EQE notebook, seen through a different observable.
    """
    W_cm = np.asarray(W_um, dtype=float) * 1e-4
    inv = 1.0 / (np.asarray(tau_bulk_us, dtype=float) * 1e-6) \
        + 2 * np.asarray(S_cm_s, dtype=float) / W_cm
    return 1e6 / inv


# --- the charge balance ----------------------------------------------------
def solve_surface_potential(Qf_cm2, doping_cm3, dopant_type='n',
                            delta_n_cm3=0.0, K_s=K_SI, T=T_ROOM, ni=NI_SI,
                            dit_midgap=1e11, dit_edge=1e13, sigma_n=1e-16,
                            sigma_p=1e-16, Eg_eV=EG_SI_EV):
    """
    Surface potential from charge neutrality - Eq. (17).

        Q_f + Q_sc(phi_s) + Q_it(phi_s) = 0

    Everything the dielectric puts on the surface has to be mirrored
    somewhere, and the only free variable is phi_s. Solved by bracketing
    and bisection over a range slightly wider than the band gap.

    Set delta_n_cm3 = 0 for the dark case and > 0 for the illuminated one;
    the difference between the two CPDs is the surface photovoltage.
    """
    def residual(phi):
        return (Qf_cm2
                + space_charge_density(phi, doping_cm3, dopant_type,
                                       delta_n_cm3, K_s, T, ni)
                + interface_charge(phi, doping_cm3, dopant_type, delta_n_cm3,
                                   dit_midgap, dit_edge, sigma_n, sigma_p,
                                   Eg_eV, T, ni))

    for lo, hi in ((-(Eg_eV + 0.3), Eg_eV + 0.3), (-3.0, 3.0)):
        if residual(lo) * residual(hi) <= 0:
            return brentq(residual, lo, hi, xtol=1e-12, rtol=1e-14)
    raise ValueError(
        f"charge balance not bracketed for Qf = {Qf_cm2:.3g} cm^-2; the "
        "dielectric charge is outside the range this model can mirror in "
        "the semiconductor")


def surface_regime(phi_s_V, dopant_type='n', flat_band_mV=10.0,
                   Eg_eV=EG_SI_EV):
    """
    Name the surface condition for a given band bending - the four cases
    tabulated in Sec. 10. On n-type material phi_s > 0 is accumulation; on
    p-type it is depletion or inversion.
    """
    phi = float(phi_s_V)
    if abs(phi) < flat_band_mV * 1e-3:
        return 'flat band'
    n_type = dopant_type.lower().startswith('n')
    majority_accumulates = (phi > 0) if n_type else (phi < 0)
    if majority_accumulates:
        return 'accumulation'
    return 'inversion' if abs(phi) > Eg_eV / 2 else 'depletion'


def cpd(phi_s_V, Qf_cm2=0.0, xc_nm=0.0, doping_cm3=5e15, dopant_type='n',
        phi_probe_eV=WORK_FUNCTION_EV['Au'], K_i=K_SIO2, chi_eV=CHI_SI_EV,
        Eg_eV=EG_SI_EV, T=T_ROOM, ni=NI_SI):
    """
    The master equation - Eq. (10).

        CPD = -( Phi_ms/q + phi_s + V_i )

    with Phi_ms = Phi_probe - Phi_semiconductor the flat-band work function
    difference in eV, phi_s the band bending, and V_i the drop across the
    dielectric from Eq. (14).

    One measurement returns one number containing all three terms, so no
    single CPD reading can separate them. Separating them is what the rest
    of the notebook is about.
    """
    phi_s_semi = work_function_semiconductor(doping_cm3, dopant_type, chi_eV,
                                             Eg_eV, T, ni)
    phi_ms = phi_probe_eV - phi_s_semi
    V_i = insulator_potential(Qf_cm2, xc_nm, K_i)
    return -(phi_ms + np.asarray(phi_s_V, dtype=float) + V_i)


# ===========================================================================
# PART 3 - light
# ===========================================================================
def surface_photovoltage(cpd_dark_V, cpd_light_V):
    """
    Surface photovoltage - Eq. (18).

        SPV = CPD_dark - CPD_light

    Because it is a difference of two CPDs taken with the same probe, the
    probe work function, the calibration and the dielectric term all cancel
    exactly. SPV is therefore reproducible to a few mV where an absolute
    work function is good to perhaps 50 meV.
    """
    return np.asarray(cpd_dark_V, float) - np.asarray(cpd_light_V, float)


def cpd_dark_and_light(Qf_cm2, delta_n_cm3, doping_cm3=5e15,
                       dopant_type='n', xc_nm=5.0, K_i=K_SIO2,
                       phi_probe_eV=WORK_FUNCTION_EV['Au'], **kwargs):
    """
    Solve the charge balance twice - dark and illuminated - and return
    (CPD_dark, CPD_light, SPV, phi_s_dark, phi_s_light).

    The convenience wrapper the notebook uses for Eqs. (17), (10) and (18)
    in sequence. Extra keyword arguments go to solve_surface_potential.
    """
    common = dict(doping_cm3=doping_cm3, dopant_type=dopant_type, **kwargs)
    phi_dark = solve_surface_potential(Qf_cm2, delta_n_cm3=0.0, **common)
    phi_light = solve_surface_potential(Qf_cm2, delta_n_cm3=delta_n_cm3,
                                        **common)
    cpd_kw = dict(Qf_cm2=Qf_cm2, xc_nm=xc_nm, doping_cm3=doping_cm3,
                  dopant_type=dopant_type, phi_probe_eV=phi_probe_eV,
                  K_i=K_i)
    c_dark = cpd(phi_dark, **cpd_kw)
    c_light = cpd(phi_light, **cpd_kw)
    return (c_dark, c_light, surface_photovoltage(c_dark, c_light),
            phi_dark, phi_light)


def charge_fluctuation_average(Qf0_cm2, sigma_q_cm2, quantity, n_sigma=3.0,
                               n_points=25):
    """
    Average a surface quantity over a Gaussian spread of dielectric charge -
    Eq. (20).

        <X> = int X(Q_f) exp(-(Q_f - Q_f0)^2 / 2 sigma_q^2) dQ_f
              / int exp(...) dQ_f

    Real dielectric films are not uniformly charged, and a millimetre-scale
    probe averages over the patchwork. The visible effect is that the sharp
    accumulation-to-inversion transition of Eq. (17) is smeared out over a
    wider range of mean charge.

    `quantity` is any callable taking a single Q_f in cm^-2.
    """
    if sigma_q_cm2 <= 0:
        return quantity(Qf0_cm2)
    q_grid = np.linspace(Qf0_cm2 - n_sigma * sigma_q_cm2,
                         Qf0_cm2 + n_sigma * sigma_q_cm2, n_points)
    w = np.exp(-0.5 * ((q_grid - Qf0_cm2) / sigma_q_cm2) ** 2)
    vals = np.array([quantity(q) for q in q_grid])
    return _trapz(vals * w, q_grid) / _trapz(w, q_grid)


# --- the SPV diffusion-length method --------------------------------------
_ALPHA_NM = np.array([700, 750, 800, 850, 900, 950, 1000, 1050, 1100])
_ALPHA_CM1 = np.array([1.9e3, 1.1e3, 6.5e2, 4.0e2, 2.4e2, 1.4e2, 7.0e1,
                       3.4e1, 1.5e1])


def alpha_silicon(wavelength_nm):
    """
    Absorption coefficient of crystalline silicon in 1/cm, log-interpolated
    from literature anchor points over 700-1100 nm - the range the SPV
    diffusion-length method uses. The same coarse model as the EQE
    notebook's helper; adequate for shapes, not for certification.
    """
    wl = np.asarray(wavelength_nm, dtype=float)
    return 10 ** np.interp(wl, _ALPHA_NM, np.log10(_ALPHA_CM1))


def spv_constant_flux(alpha_cm1, L_n_um, C2=1.0):
    """
    Small-signal SPV under constant photon flux - Eq. (23).

        1/V_P = C2 * (L_n + 1/alpha)

    Returns V_P (arbitrary units, set by C2). The useful feature is that
    1/V_P is linear in 1/alpha with an intercept on the negative 1/alpha
    axis at -L_n.
    """
    L_cm = np.asarray(L_n_um, dtype=float) * 1e-4
    return 1.0 / (C2 * (L_cm + 1.0 / np.asarray(alpha_cm1, dtype=float)))


def fit_diffusion_length(alpha_cm1, V_P, sigma=None):
    """
    The Goodman extraction - Eq. (23) rearranged.

    Fits 1/V_P against 1/alpha and returns (L_n in um, one-sigma error in
    um, (slope, intercept)). L_n is the *negative* x-intercept, so it is an
    extrapolation outside the measured range - the same slope-is-safe /
    intercept-is-fragile asymmetry as the TLM notebook, and the reason this
    method needs wavelengths whose 1/alpha straddle the expected L_n.
    """
    x = 1.0 / np.asarray(alpha_cm1, dtype=float)          # cm
    y = 1.0 / np.asarray(V_P, dtype=float)
    w = None if sigma is None else 1.0 / np.asarray(sigma, dtype=float)
    (m, b), cov = np.polyfit(x, y, 1, w=w, cov=True)
    L_cm = b / m
    var = ((1.0 / m) ** 2 * cov[1, 1] + (b / m ** 2) ** 2 * cov[0, 0]
           - 2 * (1.0 / m) * (b / m ** 2) * cov[0, 1])
    return L_cm * 1e4, np.sqrt(max(var, 0.0)) * 1e4, (m, b)


# ===========================================================================
# PART 4 - synthetic measurements, with the defects real ones have
# ===========================================================================
def synthetic_null_sweep(V_backing, cpd_true_V, gradient_V_per_V=3.0,
                         noise_mV=2.0, spacing_drift=0.0, seed=0):
    """
    An off-null sweep as a real instrument delivers it - Eq. (6).

    Adds Gaussian noise to the peak-to-peak signal and, optionally, a
    fractional drift in the mean spacing across the sweep.

    The two are not equivalent faults. A spacing that is *steady but wrong*
    changes the gradient k and leaves the zero crossing exactly where it
    was, so CPD is unaffected. A spacing that *drifts during the sweep*
    makes V_ptp a curve rather than a line, and a straight-line fit to a
    curve crosses zero in the wrong place - a systematic error in CPD that
    the residuals barely show. That is what the active spacing regulation
    on a real instrument is for.
    """
    rng = np.random.default_rng(seed)
    V_b = np.asarray(V_backing, dtype=float)
    ramp = 1.0 + spacing_drift * np.linspace(0, 1, V_b.size)
    signal = gradient_V_per_V * ramp * (V_b + cpd_true_V)
    return signal + rng.normal(0.0, noise_mV * 1e-3, V_b.size)


def synthetic_probe_drift(minutes, cpd_true_V, drift_mV_per_hour=25.0,
                          noise_mV=1.0, seed=1):
    """
    Repeated CPD readings on one unchanging sample over a session, with the
    slow drift of a contaminating probe superimposed. The reason absolute
    work functions are bracketed by reference measurements.
    """
    rng = np.random.default_rng(seed)
    t = np.asarray(minutes, dtype=float)
    return (cpd_true_V + drift_mV_per_hour * 1e-3 * t / 60.0
            + rng.normal(0.0, noise_mV * 1e-3, t.size))


def synthetic_spv_intensity(delta_n_cm3, Qf_cm2, doping_cm3=5e15,
                            dopant_type='n', noise_mV=1.5, seed=2, **kwargs):
    """
    SPV against injection level, solved properly at each point and then
    given measurement noise. Sweep delta_n over several decades to see the
    approach to flat band - and to see what happens if you stop too early.
    """
    rng = np.random.default_rng(seed)
    dn = np.atleast_1d(np.asarray(delta_n_cm3, dtype=float))
    out = np.array([cpd_dark_and_light(Qf_cm2, d, doping_cm3, dopant_type,
                                       **kwargs)[2] for d in dn])
    return out + rng.normal(0.0, noise_mV * 1e-3, out.size)


def synthetic_spv_wavelength_scan(wavelength_nm, L_n_um, C2=1.0,
                                  noise_frac=0.02, seed=3):
    """
    A constant-photon-flux SPV wavelength scan for the diffusion-length
    method - Eq. (23) with multiplicative noise. Returns V_P.
    """
    rng = np.random.default_rng(seed)
    wl = np.asarray(wavelength_nm, dtype=float)
    V = spv_constant_flux(alpha_silicon(wl), L_n_um, C2)
    return V * (1.0 + rng.normal(0.0, noise_frac, wl.size))
