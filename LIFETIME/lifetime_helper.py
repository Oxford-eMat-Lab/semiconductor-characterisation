"""
lifetime_helper.py
------------------
Helper functions for the effective-lifetime teaching notebook,
`lifetime_analysis.ipynb`. Equation numbers in the docstrings, e.g.
"Eq. (5)", refer to the numbered equations in that notebook.

Contents
--------
Part I  (Eqs. 1-44)  Physics of the quantity being measured: excess
        carriers, the generation-recombination balance, SRH / radiative /
        Auger / surface recombination, how they combine into an effective
        lifetime, and the implied voltage that lifetime allows.
Part II (Eqs. 45-57) The measurement: how a Sinton-type photoconductance
        instrument turns a coil voltage and a reference-cell voltage into
        excess carrier density and generation, and how lifetime is
        extracted, checked and reported.
Embedded data (end of file): a coil-calibration data set, two example
        instrument exports and the instrument schematic, so the notebook
        needs no other file.

The Sinton instrument senses photoconductance. Lifetime is inferred from
a carrier balance; it is never measured directly.

These functions are written for teaching, not for precision silicon
metrology:
  - the examples assume uniform, quasi-neutral, nondegenerate silicon,
    full dopant ionisation and no bandgap narrowing,
  - the combined teaching-model temperature range is 275-363 K,
  - bulk SRH uses a generic independent trap with specified capture times,
  - the basic surface model uses specified effective velocities; the
    interface-state example takes a prescribed surface potential,
  - optical recycling, Schenk narrowing and incomplete ionisation are
    identified extensions, not silently mixed into every plot,
  - the measurement analysis is not a replacement for the instrument
    manufacturer's calibration or analysis software.

UNIT AND SIGN CONVENTIONS
-------------------------
Carrier/dopant densities: cm^-3. Volume recombination: cm^-3 s^-1.
Surface flux: cm^-2 s^-1. Lifetime: s internally, converted to us on plots.
Temperature: K. Energies: eV. Potential and implied voltage: V.
Sample.W_um is in um; Sample.W_cm converts thickness to cm internally.
Capture cross sections: cm^2. Interface-state density: cm^-2 eV^-1.
Doping is a positive magnitude with an explicit kind, 'p' or 'n'.
Surface potential follows n_s = n_b exp(psi_s/Vt).
Part II: lengths in cm, except coil geometry (lambda, lift-off, calibration
thickness), which is in mm and says so; sheet conductance in S; time in s.
"""

from dataclasses import dataclass, replace
import base64
import hashlib
import io
import zlib

import numpy as np
import pandas as pd
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Physical constants (cm-based semiconductor quantities; energies in eV)
# ---------------------------------------------------------------------------
Q = 1.602176634e-19  # elementary charge, C
KB_EV = 8.617333262145e-5  # Boltzmann constant, eV/K
KB_J = 1.380649e-23  # Boltzmann constant, J/K
M0 = 9.1093837015e-31  # electron rest mass, kg
EPS_SI = 11.7 * 8.8541878128e-14  # F/cm
# NumPy 2 renamed trapz; support both NumPy generations.
TRAPZ = getattr(np, "trapezoid", None) or np.trapz


# ---------------------------------------------------------------------------
# Input validation and sample parameters
# ---------------------------------------------------------------------------
def _check(value, name, positive=False):
    """Validate finite, nonnegative (or positive) physical inputs."""
    a = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(a)) or np.any(a <= 0 if positive else a < 0):
        raise ValueError(
            f'{name} must be finite and {"positive" if positive else "nonnegative"}.'
        )
    return a


def _temperature(T, lo=275.0, hi=363.0):
    """Enforce the stated common teaching-model temperature interval, in K."""
    T = _check(T, "temperature", True)
    if np.any((T < lo) | (T > hi)):
        raise ValueError(f"Temperature must lie within {lo:g}–{hi:g} K for this model.")
    return T


@dataclass(frozen=True)
class Sample:
    """Uniform wafer inputs for Eqs. (1), (35).

    Doping is a magnitude in cm^-3; kind is 'p' or 'n'. W_um is wafer thickness
    in micrometres. S_front/S_rear are effective velocities in cm/s. Capture
    times are in seconds, trap_energy is E_t-E_i in eV. f_reabs is the assumed
    band-to-band photon recycling probability, not the Sinton optical factor.
    """

    doping: float = 1e16
    kind: str = "p"
    T: float = 300.0
    W_um: float = 180.0
    S_front: float = 5.0
    S_rear: float = 5.0
    tau_n0: float = 2e-3
    tau_p0: float = 8e-3
    trap_energy: float = 0.0
    f_reabs: float = 0.0

    def __post_init__(self):
        """Check sample parameters against the teaching model domain."""
        _temperature(self.T)
        for key in ("doping", "S_front", "S_rear"):
            _check(getattr(self, key), key)
        for key in ("W_um", "tau_n0", "tau_p0"):
            _check(getattr(self, key), key, True)
        if self.kind not in ("p", "n"):
            raise ValueError("kind must be 'p' or 'n'.")
        if not np.isfinite(self.trap_energy) or abs(self.trap_energy) > 0.45:
            raise ValueError("Use a trap within 0.45 eV of E_i in these examples.")
        if not 0 <= self.f_reabs < 1:
            raise ValueError("f_reabs must lie in [0, 1).")

    @property
    def W_cm(self):
        """Wafer thickness in cm, used in Eq. (35)."""
        return self.W_um * 1e-4


# ---------------------------------------------------------------------------
# 1. Bandgap, intrinsic density and injected carrier populations
# ---------------------------------------------------------------------------
def bandgap_passler(T):
    """Unperturbed Si gap in eV, Pässler (2002), Eq. (26).

    Temperature in K; restricted here to 275–363 K, the range used by
    every model in this module.
    """
    T = _temperature(T)
    alpha, delta, theta = 0.323e-3, 0.51, 446.0
    x = 2 * T / theta
    polynomial = (
        1
        + np.pi**2 / (3 * (1 + delta**2)) * x**2
        + (0.75 * delta**2 - 0.25) * x**3
        + 8 / 3 * x**4
        + x**6
    )
    return 1.17 - alpha * theta * (
        (1 - 3 * delta**2) / np.expm1(theta / T)
        + 1.5 * delta**2 * (polynomial ** (1 / 6) - 1)
    )


def intrinsic_density(T):
    """Unperturbed ni in cm^-3, Couderc et al. (2014), Eq. (27).

    Uses the Pässler gap. This is the 1.541e15*T**1.712 form of Couderc et al.
    """
    T = _temperature(T)
    return 1.541e15 * T**1.712 * np.exp(-bandgap_passler(T) / (2 * KB_EV * T))


def carrier_populations(delta_n, doping=1e16, kind="p", T=300.0, ni=None):
    """Return n, p, n0, p0, ni (cm^-3), Eqs. (1), (3).

    Nondegenerate, fully ionised, uncompensated wafer; uniform paired excess.
    ni can be supplied explicitly for a controlled statistics comparison.
    The minority population uses mass action to avoid subtractive cancellation.
    """
    dn = _check(delta_n, "excess density")
    N = _check(doping, "doping")
    _temperature(T)
    if kind not in ("p", "n"):
        raise ValueError("kind must be 'p' or 'n'.")
    ni = intrinsic_density(T) if ni is None else _check(ni, "ni", True)
    major = 0.5 * (N + np.hypot(N, 2 * ni))
    minor = ni**2 / major
    n0, p0 = (minor, major) if kind == "p" else (major, minor)
    return n0 + dn, p0 + dn, n0, p0, ni


# ---------------------------------------------------------------------------
# 2. Lifetime from net loss and independent rate addition
# ---------------------------------------------------------------------------
def lifetime_from_rate(delta_n, rate):
    """Lifetime Δn/U in seconds, Eq. (5).

    Zero net rate at positive injection gives infinity. At thermal equilibrium
    the ratio is undefined and returns NaN. Negative rates are rejected.
    """
    dn, rate = np.broadcast_arrays(
        _check(delta_n, "excess density"), _check(rate, "rate")
    )
    out = np.full(dn.shape, np.nan)
    np.divide(dn, rate, out=out, where=rate > 0)
    return np.where((rate == 0) & (dn > 0), np.inf, out)


def combine_lifetimes(*lifetimes):
    """Harmonic combination of lifetimes in seconds, Eq. (36).

    Infinite lifetime represents an absent channel; finite values must be
    positive. All lifetimes must refer to the same carrier population.
    """
    if not lifetimes:
        raise ValueError("Supply at least one lifetime.")
    values = np.broadcast_arrays(*[np.asarray(x, dtype=float) for x in lifetimes])
    if any(np.any(np.isnan(x) | (x <= 0)) for x in values):
        raise ValueError("Lifetimes must be positive or +infinity.")
    inv = sum(1 / x for x in values)
    return np.divide(1.0, inv, out=np.full(inv.shape, np.inf), where=inv > 0)


# ---------------------------------------------------------------------------
# 3. Bulk recombination: SRH, radiative and Black Auger models
# ---------------------------------------------------------------------------
def capture_times(N_t, sigma_n, sigma_p, T=300.0, mass_n=0.26, mass_p=0.39):
    """Electron/hole capture times in s, Eq. (9).

    N_t in cm^-3, cross sections in cm^2, effective masses in units of m0.
    Uses mean Maxwell speed sqrt(8 kT/(pi m)); distinguish this convention
    from an RMS thermal speed when importing capture cross sections.
    """
    T = _temperature(T)
    Nt, sn, sp = [
        _check(x, name, True)
        for x, name in [(N_t, "N_t"), (sigma_n, "sigma_n"), (sigma_p, "sigma_p")]
    ]
    mn, mp = _check(mass_n, "mass_n", True), _check(mass_p, "mass_p", True)
    vn, vp = (
        np.sqrt(8 * KB_J * T / (np.pi * mn * M0)) * 100,
        np.sqrt(8 * KB_J * T / (np.pi * mp * M0)) * 100,
    )
    return 1 / (Nt * sn * vn), 1 / (Nt * sp * vp)


def srh_rate(delta_n, sample):
    """Single-level SRH net rate in cm^-3 s^-1, Eqs. (8), (10).

    delta_n : excess carrier density, cm^-3 (scalar or array).
    sample : Sample; tau_n0 and tau_p0 are capture times in seconds.
    Returns the net bulk SRH rate, cm^-3 s^-1, shaped like delta_n.

        U_SRH = (np - ni^2) / [tau_p0 (n+n1) + tau_n0 (p+p1)]

    Trap energy is E_t-E_i, not E_t-E_v. The excess product is evaluated as
    Δn(n0+p0+Δn), exactly equal to np-ni^2 in this model but stable near zero.
    """
    dn = _check(delta_n, "excess density")
    n, p, n0, p0, ni = carrier_populations(dn, sample.doping, sample.kind, sample.T)
    n1 = ni * np.exp(sample.trap_energy / (KB_EV * sample.T))
    p1 = ni * np.exp(-sample.trap_energy / (KB_EV * sample.T))
    return dn * (n0 + p0 + dn) / (sample.tau_p0 * (n + n1) + sample.tau_n0 * (p + p1))


def radiative_coefficient(T):
    """B_low in cm^3/s, Nguyen et al. (2014) Eq. (5); notebook Eq. (12).

    The source's B coefficient is tied to its intrinsic-density convention.
    The 275–363 K restriction is the common range used in this notebook.
    """
    T = _temperature(T)
    return 10 ** np.polynomial.polynomial.polyval(
        T, [-9.65614, -8.05258e-2, 6.02695e-4, -2.29844e-6, 4.31934e-9, -3.16154e-12]
    )


def radiative_rate(delta_n, sample, B=None):
    """Net radiative loss in cm^-3 s^-1, Eqs. (11), (14).

    Constant B (if given) or temperature-dependent B_low. The main examples
    omit BGN/B_rel and use f_reabs as an explicit controlled input, not a
    calibrated optical result. Recycled pairs are excluded from net loss.
    """
    dn = _check(delta_n, "excess density")
    _, _, n0, p0, _ = carrier_populations(dn, sample.doping, sample.kind, sample.T)
    B = radiative_coefficient(sample.T) if B is None else _check(B, "B")
    return B * (1 - sample.f_reabs) * dn * (n0 + p0 + dn)


def auger_coefficients(n, p, T=300.0, model="black"):
    """Effective C_n,C_p in cm^6/s; Eqs. (18), (19).

    'black': Black et al. (2026) parameter set.
    'constant': screened 300 K coefficients 3.41e-31 and 1.17e-31, used only
    to expose the effect of omitting Coulomb enhancement (not a competing
    complete intrinsic model). The Black function alone allows 117–463 K.
    """
    n, p = _check(n, "n"), _check(p, "p")
    T = _temperature(T, 117.0, 463.0)
    if model == "constant":
        return np.zeros_like(n + p) + 3.41e-31, np.zeros_like(n + p) + 1.17e-31
    if model != "black":
        raise ValueError("model must be 'black' or 'constant'.")
    cn = 1.02e-31 * np.exp(-2.8e-3 / (KB_EV * T)) + 9.71e-32 / np.tanh(
        21.2e-3 / (2 * KB_EV * T)
    )
    cp = 4.51e-32 / np.tanh(21.0e-3 / (2 * KB_EV * T))
    temp_enh = 1 + 880.9 * np.exp(-((T / 17.51) ** 0.6434))
    density = 1 + ((n + p) / (4e17 * (T / 300) ** 4)) ** 2
    return cn * (1 + (1.629 * temp_enh - 1) / density), cp * (
        1 + (1.815 * temp_enh - 1) / density
    )


def auger_rate(delta_n, sample, model="black"):
    """Net Auger rate in cm^-3 s^-1, Eq. (17).

    Uses the separate cubic equilibrium subtractions in the Black model.
    Expanded polynomials avoid cancellation at low Δn. Coefficients are
    evaluated at illuminated populations, as in that parameterisation.
    """
    dn = _check(delta_n, "excess density")
    n, p, n0, p0, _ = carrier_populations(dn, sample.doping, sample.kind, sample.T)
    cn, cp = auger_coefficients(n, p, sample.T, model)
    eeh = dn * (n0 * n0 + 2 * n0 * p0) + dn**2 * (2 * n0 + p0) + dn**3
    ehh = dn * (p0 * p0 + 2 * n0 * p0) + dn**2 * (2 * p0 + n0) + dn**3
    return cn * eeh + cp * ehh


def component_lifetimes(delta_n, sample=Sample(), auger_model="black"):
    """Component, bulk, intrinsic and effective lifetimes in s, Eqs. (36), (38), (35).

    delta_n : positive excess carrier density, cm^-3 (scalar or array).
    sample : Sample containing doping, temperature, thickness and trap inputs.
    auger_model : 'black' or the explicitly simplified 'constant' comparison.

    Returns a dict of lifetimes in seconds, with keys srh, radiative, auger,
    surface, intrinsic, bulk and effective, all shaped like delta_n.
    The rate sum gives tau_eff = 1 / sum(1 / tau_j), not an arithmetic mean.

    Uniform carrier density; equivalent density used for all channels. No
    mobility, calibration, trapping or carrier-profile model is included.
    """
    dn = _check(delta_n, "excess density", True)
    srh = lifetime_from_rate(dn, srh_rate(dn, sample))
    rad = lifetime_from_rate(dn, radiative_rate(dn, sample))
    aug = lifetime_from_rate(dn, auger_rate(dn, sample, auger_model))
    surf = np.full_like(
        dn,
        (
            np.inf
            if sample.S_front + sample.S_rear == 0
            else sample.W_cm / (sample.S_front + sample.S_rear)
        ),
    )
    intrinsic = combine_lifetimes(rad, aug)
    bulk = combine_lifetimes(srh, intrinsic)
    return dict(
        srh=srh,
        radiative=rad,
        auger=aug,
        surface=surf,
        intrinsic=intrinsic,
        bulk=bulk,
        effective=combine_lifetimes(bulk, surf),
    )


# ---------------------------------------------------------------------------
# 4. Carrier storage during illumination and decay
# ---------------------------------------------------------------------------
def rectangular_generation_response(time_s, G, tau_s, pulse_s):
    """Analytic constant-lifetime population response in cm^-3, Eq. (7).

    G is a uniform pair-generation rate in cm^-3 s^-1. Illumination starts
    at t=0 and ends at pulse_s; this is a balance illustration, not a Sinton
    flash simulation or mode-analysis implementation.
    """
    t = _check(time_s, "time")
    G, tau, on = (
        _check(G, "G"),
        _check(tau_s, "tau", True),
        _check(pulse_s, "pulse", True),
    )
    return (
        G
        * tau
        * (-np.expm1(-np.minimum(t, on) / tau))
        * np.exp(-np.maximum(t - on, 0) / tau)
    )


# ---------------------------------------------------------------------------
# 5. Implied voltage and the generation-defined operating point
# ---------------------------------------------------------------------------
def implied_voltage(delta_n, sample=Sample(), ni=None):
    """Quasi-Fermi splitting/q in V, Eq. (43).

    Nondegenerate and spatially uniform; no contact/series-resistance losses.
    Uses log1p for zero injection. A value at arbitrary injection is not
    automatically one-sun implied Voc; solve the generation balance first.
    """
    dn = _check(delta_n, "excess density")
    _, _, n0, p0, ni = carrier_populations(dn, sample.doping, sample.kind, sample.T, ni)
    return KB_EV * sample.T * np.log1p(dn * (n0 + p0 + dn) / ni**2)


def operating_point(sample=Sample(), Jgen_mA_cm2=40.0):
    """Solve uniform open-circuit generation balance, Eq. (44).

    Jgen is the absorbed pair-generation equivalent current density, not an
    independently guaranteed cell Jsc. Returns Δn (cm^-3), lifetime (s), and
    implied voltage (V). No extrapolation outside Δn=1e8–1e18 cm^-3.
    """
    G = float(_check(Jgen_mA_cm2, "Jgen", True)) * 1e-3 / (Q * sample.W_cm)

    def residual(log_dn):
        """Net loss minus generation at the trial log10 excess density."""
        dn = 10**log_dn
        return dn / component_lifetimes(dn, sample)["effective"] - G

    dn = 10 ** brentq(residual, 8, 18)
    return dict(
        delta_n=dn,
        tau=float(component_lifetimes(dn, sample)["effective"]),
        voltage=float(implied_voltage(dn, sample)),
        G=G,
    )


# ---------------------------------------------------------------------------
# 6. Radiative density correction, recycling and absorption
# ---------------------------------------------------------------------------
def b_relative(n, p, T=300.0):
    """Dimensionless Altermatt B/B_low extension, Eq. (13).

    Uses total n+p with doubled b1,b3 scales.
    This illustrates a carrier-density correction and is not silently folded
    into the unperturbed main examples. Match its ni/BGN convention before
    quantitative use (Fell et al., 2021).
    """
    T = _temperature(T)
    total = _check(n, "n") + _check(p, "p")
    bmin = 0.20 - 0.20 / (1 + (T / 320) ** 2.50)
    b1 = 2 * (1.5e18 + (1e7 - 1.5e18) / (1 + (T / 550) ** 3))
    b3 = 2 * (4e18 + (1e9 - 4e18) / (1 + (T / 365) ** 3.54))
    return bmin + (1 - bmin) / (1 + (total / b1) ** 0.54 + (total / b3) ** 1.25)


def photon_recycling_loss(wavelength_nm, emission_weight, alpha_bb, alpha_fca, p_abs):
    """Return net radiative loss fraction 1-f_PR, Eq. (14).

    Wavelength in nm, absorption coefficients in cm^-1. p_abs is the total
    absorption probability for the specified geometry, given by an optical
    calculation; this helper does not infer planar/diffuse light trapping.
    The weighting must be a density per wavelength, not per photon energy.
    """
    wl = _check(wavelength_nm, "wavelength", True)
    arrays = [
        np.asarray(x, float) for x in (emission_weight, alpha_bb, alpha_fca, p_abs)
    ]
    if (
        wl.ndim != 1
        or wl.size < 2
        or np.any(np.diff(wl) <= 0)
        or any(x.shape != wl.shape for x in arrays)
    ):
        raise ValueError(
            "Provide equal-length one-dimensional arrays on an increasing grid."
        )
    em, bb, fca, pa = arrays
    for x in arrays:
        _check(x, "optical input")
    if np.any(pa > 1) or TRAPZ(em, wl) <= 0:
        raise ValueError(
            "Absorption probability must be <=1 and emission weight nonzero."
        )
    bb_fraction = np.divide(bb, bb + fca, out=np.zeros_like(bb), where=bb + fca > 0)
    return float(1 - TRAPZ(em * pa * bb_fraction, wl) / TRAPZ(em, wl))


def absorption_temperature(alpha_295, c_T_per_K, T):
    """Absorption in cm^-1 from local Schinke-table inputs, Eq. (15).

    c_T is the relative coefficient in K^-1; the dimensionless power exponent
    is 295*c_T, not c_T alone. No optical table is bundled or extrapolated.
    """
    T = _temperature(T)
    return _check(alpha_295, "absorption") * (T / 295) ** (
        295 * np.asarray(c_T_per_K, float)
    )


def free_carrier_absorption(wavelength_nm, density, dopant="P"):
    """Baker-Finch FCA in cm^-1, Eq. (16).

    Enforces the paper's 1000–1500 nm and 1e18–3e20 cm^-3 range; density is
    the doped-region carrier/dopant density for the chosen P/B calibration.
    This is not a validated arbitrary injected-plasma absorption model.
    """
    wl, N = _check(wavelength_nm, "wavelength", True), _check(density, "density", True)
    if np.any((wl < 1000) | (wl > 1500)) or np.any((N < 1e18) | (N > 3e20)):
        raise ValueError("Baker-Finch range: 1000–1500 nm and 1e18–3e20 cm^-3.")
    if dopant not in ("P", "B"):
        raise ValueError("dopant must be 'P' or 'B'.")
    C, power = (1.68e-6, 2.88) if dopant == "P" else (1.82e-9, 2.18)
    return C * N * (wl * 1e-7) ** power


# ---------------------------------------------------------------------------
# 7. Carrier-statistics extensions: narrowing and ionisation
# ---------------------------------------------------------------------------
def schenk_gap_narrowing(n, p, N_D, N_A, T=300.0):
    """Positive Si gap reduction in eV, Eqs. (28), (31).

    Schenk (1998) approximation. Input carrier
    and ionised dopant densities are cm^-3. Evaluates one specified state;
    it does not iterate Poisson/statistics or establish degenerate statistics.
    """
    T = _temperature(T)
    n, p, nd, na = np.broadcast_arrays(
        *[_check(x, "density") * (3.719e-7) ** 3 for x in (n, p, N_D, N_A)]
    )
    F = KB_EV * T / 0.01655
    nx = n + p
    ni = nd + na
    px = 0.5187 * n + 0.4813 * p
    pi = 0.5187 * nd + 0.4813 * na
    ui = (ni / F) ** 2
    shift = np.zeros_like(nx)
    for carrier, g, a, b, c, d, h, j, k, z in [
        (n, 12.0, 0.5187, 8.0, 1.3346, 0.893, 3.91, 2.8585, 0.012, 0.75),
        (p, 4.0, 0.4813, 1.0, 1.2365, 1.153, 4.2, 2.9307, 0.19, 0.25),
    ]:
        xc_num = (
            (4 * np.pi) ** 3
            * nx**2
            * (
                (48 * carrier / (np.pi * g)) ** (1 / 3)
                + c * np.log1p(d * px ** (7 / 30))
            )
            + 8 * np.pi * a / g * carrier * F**2
            + np.sqrt(8 * np.pi * nx) * F**2.5
        )
        xc_den = (
            (4 * np.pi) ** 3 * nx**2 + F**3 + b * np.sqrt(nx) * F**2 + 40 * nx**1.5 * F
        )
        ion_den = np.sqrt(F * ni / (2 * np.pi)) * (
            1 + h * np.log1p(np.sqrt(ni) / F)
        ) + j * ui * pi**0.75 * (1 + k * pi**z)
        ion = np.divide(
            ni * (1 + ui), ion_den, out=np.zeros_like(nx), where=ion_den > 0
        )
        shift += 0.01655 * (xc_num / xc_den + ion)
    return shift


def incomplete_ionisation(
    N_dop, T=300.0, kind="p", binding_eV=0.045, Nc=2.86e19, Nv=3.10e19
):
    """Equilibrium majority/minority carriers in cm^-3, Eq. (32).

    Isolated shallow dopants, no BGN, no compensation, Boltzmann statistics.
    Nc,Nv are explicit 300 K teaching inputs scaled as T^1.5. This avoids
    silently mixing an independent fitted DOS with the chosen ni convention.
    The majority density is solved by neutrality; minority follows mass action.
    """
    T = float(_temperature(T))
    N = float(_check(N_dop, "doping", True))
    ni = float(intrinsic_density(T))
    if kind not in ("p", "n"):
        raise ValueError("kind must be 'p' or 'n'.")
    dos = (Nv if kind == "p" else Nc) * (T / 300) ** 1.5
    factor = (4 if kind == "p" else 2) * np.exp(binding_eV / (KB_EV * T)) / dos

    def f(x):
        """Charge-neutrality residual for a trial majority density."""
        return x - ni**2 / x - N / (1 + factor * x)

    majority = brentq(f, ni, N + ni, xtol=1e-3)
    return majority, ni**2 / majority


# ---------------------------------------------------------------------------
# 8. Interface states and surface recombination
# ---------------------------------------------------------------------------
def interface_distribution(
    E, peak=1e11, centre=0.0, width=0.20, kind="gaussian", endpoint=1e8
):
    """Dit(E) in cm^-2 eV^-1, Eq. (20).

    E and centre are relative to E_i in eV. width is sigma for a Gaussian,
    full width for a tophat, or signed energy interval for a tail. Tail
    support must be selected by the caller; the function defines its shape.
    """
    E = np.asarray(E, float)
    peak = float(_check(peak, "peak", True))
    if not np.isfinite(width) or width == 0:
        raise ValueError("width must be finite and nonzero.")
    if kind == "gaussian":
        return peak * np.exp(-0.5 * ((E - centre) / abs(width)) ** 2)
    if kind == "tophat":
        return peak * (np.abs(E - centre) < abs(width) / 2)
    if kind == "tail":
        return peak * np.exp(
            np.log(float(_check(endpoint, "endpoint", True)) / peak)
            * (E - centre)
            / width
        )
    raise ValueError("kind must be gaussian, tophat or tail.")


def interface_recombination(
    delta_n, psi_s, E, Dit, sigma_n=1e-15, sigma_p=1e-16, sample=Sample()
):
    """Surface SRH flux in cm^-2 s^-1, Eqs. (21), (22).

    Prescribed surface potential psi_s in V; E relative to E_i (eV), Dit per
    cm²/eV. No fitting or extraction of Q_f is implied. Capture coefficients
    are sigma*v_th (cm³/s), before multiplying by Dit in the integral.
    """
    dn = float(_check(delta_n, "excess density", True))
    E = np.asarray(E, float)
    Dit = _check(Dit, "Dit")
    if E.ndim != 1 or E.size < 2 or E.shape != Dit.shape or np.any(np.diff(E) <= 0):
        raise ValueError("E and Dit must share an increasing one-dimensional grid.")
    n, p, n0, p0, ni = carrier_populations(dn, sample.doping, sample.kind, sample.T)
    vt = KB_EV * sample.T
    psi = np.asarray(psi_s, float)
    if np.any(~np.isfinite(psi)) or np.any(np.abs(psi) > 0.30):
        raise ValueError(
            "Prescribed-potential illustration is restricted to |psi_s|<=0.30 V."
        )
    ns = n * np.exp(psi[..., None] / vt)
    ps = p * np.exp(-psi[..., None] / vt)
    vn = np.sqrt(8 * KB_J * sample.T / (np.pi * 0.26 * M0)) * 100
    vp = np.sqrt(8 * KB_J * sample.T / (np.pi * 0.39 * M0)) * 100
    cn = _check(sigma_n, "sigma_n", True) * vn
    cp = _check(sigma_p, "sigma_p", True) * vp
    n1 = ni * np.exp(E / vt)
    p1 = ni * np.exp(-E / vt)
    rate = Dit * cn * cp * dn * (n0 + p0 + dn) / (cn * (ns + n1) + cp * (ps + p1))
    return TRAPZ(rate, E, axis=-1)


def charge_average(charge_cm2, rate, mean_cm2, sigma_cm2):
    """Normalised finite-grid charge-distribution average, Eq. (24).

    Charge is elementary charges/cm²; rate is a previously calculated surface
    flux for each charge. The surface electrostatic solution is an input.
    """
    x = np.asarray(charge_cm2, float)
    r = _check(rate, "rate")
    s = float(_check(sigma_cm2, "sigma", True))
    if x.ndim != 1 or x.size < 2 or x.shape != r.shape or np.any(np.diff(x) <= 0):
        raise ValueError("Charge and rate need an increasing common grid.")
    weight = np.exp(-0.5 * ((x - mean_cm2) / s) ** 2)
    norm = TRAPZ(weight, x)
    if norm <= 0:
        raise ValueError("Charge grid misses the distribution.")
    return float(TRAPZ(weight * r, x) / norm)


def surface_J0(delta_n, U_s, sample=Sample()):
    """Surface saturation-current parameter in A/cm², Eq. (25).

    Definition at the stated injection; may depend on injection. It is not
    automatically a constant diode parameter. Nondegenerate surface/bulk
    quasi-Fermi levels are assumed equal. Zero injection is undefined.
    """
    dn = _check(delta_n, "excess density", True)
    _, _, n0, p0, ni = carrier_populations(dn, sample.doping, sample.kind, sample.T)
    return Q * _check(U_s, "surface flux") * ni**2 / (dn * (n0 + p0 + dn))


# Plotting functions: keep notebook cells short, with interpretation in Markdown.
# ---------------------------------------------------------------------------
# 9. Teaching plots and sample-condition annotations
# ---------------------------------------------------------------------------
def _figure(ncols=1, width=None, *, nrows=1, height=3.9, sharey=False):
    """Create readable Matplotlib axes using the EQE/TLM visual conventions."""
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "lines.linewidth": 1.8,
            "legend.framealpha": 0.95,
            "grid.linewidth": 0.6,
            "axes.spines.right": False,
            "figure.dpi": 120,
            "savefig.dpi": 300,
        }
    )
    return plt.subplots(
        nrows,
        ncols,
        figsize=(width or (6.8 if ncols == 1 else 11.6), height),
        sharey=sharey,
        layout="constrained",
    )


def _scientific(value):
    """Format a nonnegative plot value as a compact mathtext number."""
    if value == 0:
        return "0"
    coefficient, exponent = f"{value:.6e}".split("e")
    factor = float(coefficient)
    exponent = int(exponent)
    return (
        rf"10^{{{exponent}}}" if factor == 1 else rf"{factor:g}\times10^{{{exponent}}}"
    )


def _sample_note(
    ax,
    kind,
    doping,
    temperature,
    *,
    xy=(0.04, 0.05),
    extra=None,
    va="bottom",
    ha="left",
    fontsize=9,
):
    """Write actual sample conditions inside an axes in a compact white box.

    Doping is a positive density in cm^-3 or a description for a sweep.
    Temperature is in K or a description of its axis/colour-bar encoding.
    This is a plotting annotation; it does not change any model inputs.
    """
    symbol = "N_A" if kind == "p" else "N_D"
    density = (
        doping
        if isinstance(doping, str)
        else rf"${symbol}={_scientific(doping)}\ \mathrm{{cm}}^{{-3}}$"
    )
    temp = temperature if isinstance(temperature, str) else f"T = {temperature:g} K"
    label = f"{kind}-type Si\n{density}\n{temp}"
    if extra:
        label += "\n" + extra
    return ax.text(
        *xy,
        label,
        transform=ax.transAxes,
        va=va,
        ha=ha,
        multialignment="left",
        fontsize=fontsize,
        linespacing=1.35,
        zorder=10,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.92, pad=3),
    )


def plot_populations(doping=1e16):
    """Carrier populations and injection crossover, Eqs. (1), (3)."""
    fig, axs = _figure(2)
    dn = np.logspace(11, 18, 250)
    for ax, kind in zip(axs, ("p", "n")):
        n, p, *_ = carrier_populations(dn, doping, kind)
        ax.loglog(dn, n, label="electrons, n")
        ax.loglog(dn, p, "--", label="holes, p")
        ax.axvline(doping, color=".4", ls=":", label="Δn = doping")
        ax.set(
            xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
            ylabel="Carrier density (cm$^{-3}$)",
            title=f"{kind}-type, 300 K",
        )
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(alpha=0.2)
        _sample_note(ax, kind, doping, 300.0, xy=(0.96, 0.05), ha="right")
    return fig


def plot_population_balance():
    """Rise and fall under a rectangular generation pulse, Eq. (7)."""
    fig, ax = _figure()
    t = np.linspace(0, 0.012, 601)
    G = 1e18
    for tau in (0.0002, 0.001, 0.003):
        dn = rectangular_generation_response(t, G, tau, 0.006)
        ax.plot(t * 1e3, dn / 1e15, label=f"τ = {tau*1e3:g} ms")
    ax.axvspan(0, 6, color="gold", alpha=0.12, label="illumination on")
    ax.set(
        xlabel="Time (ms)",
        ylabel="Excess density (10$^{15}$ cm$^{-3}$)",
        title="Generation fills the population; recombination empties it",
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.text(
        0.97,
        0.55,
        "Constant-lifetime illustration\nDoping / type / T not specified",
        transform=ax.transAxes,
        ha="right",
        multialignment="left",
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.92, pad=3),
    )
    return fig


def plot_channels(sample=Sample()):
    """Bulk channel lifetimes and their sum, Eqs. (36), (38)."""
    fig, ax = _figure()
    dn = np.logspace(12, 17.5, 220)
    data = component_lifetimes(dn, sample)
    for key, lab in [
        ("srh", "SRH"),
        ("radiative", "radiative, no recycling"),
        ("auger", "Auger (Black)"),
        ("intrinsic", "intrinsic"),
        ("bulk", "total bulk"),
    ]:
        ax.loglog(
            dn,
            data[key] * 1e6,
            label=lab,
            lw=2.3 if key == "bulk" else 1.5,
            ls="--" if key == "intrinsic" else "-",
        )
    ax.set(
        xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
        ylabel="Lifetime (µs)",
        title=f"{sample.kind}-type, {sample.doping:.0e} cm$^{{-3}}$, {sample.T:g} K",
    )
    ax.legend(fontsize=8.5, ncol=2, loc="lower left", bbox_to_anchor=(0.02, 0.28))
    ax.grid(alpha=0.2)
    _sample_note(ax, sample.kind, sample.doping, sample.T)
    return fig


def plot_srh_injection():
    """SRH injection dependence without changing trap concentration, Eq. (8)."""
    fig, ax = _figure()
    dn = np.logspace(12, 17, 220)
    for tp in (0.0002, 0.002, 0.02):
        s = Sample(tau_n0=0.002, tau_p0=tp)
        ax.loglog(
            dn,
            lifetime_from_rate(dn, srh_rate(dn, s)) * 1e6,
            label=rf"$\tau_{{p0}}={tp*1e3:g}$ ms",
        )
    ax.set(
        xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
        ylabel="SRH lifetime (µs)",
        title=r"Fixed $\tau_{n0}=2$ ms, varying hole capture",
    )
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.2)
    ax.set_ylim(bottom=ax.get_ylim()[0] / 2.2)
    _sample_note(ax, s.kind, s.doping, s.T, xy=(0.96, 0.05), ha="right")
    return fig


def plot_surfaces(sample=Sample()):
    """Passivation and thickness effects in Eq. (35)."""
    fig, axs = _figure(2)
    dn = np.logspace(12, 17, 220)
    for S in (1, 10, 100):
        d = component_lifetimes(dn, replace(sample, S_front=S, S_rear=S))
        axs[0].loglog(dn, d["effective"] * 1e6, label=f"S = {S:g} cm/s per face")
    for W in (100, 180, 300):
        d = component_lifetimes(dn, replace(sample, W_um=W))
        axs[1].loglog(dn, d["effective"] * 1e6, label=f"W = {W:g} µm")
    for ax, title in zip(
        axs,
        ("Same bulk and thickness; vary surfaces", "Same bulk and S; vary thickness"),
    ):
        ax.set(
            xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
            ylabel="Effective lifetime (µs)",
            title=title,
        )
        ax.set_ylim(bottom=ax.get_ylim()[0] * 0.65, top=ax.get_ylim()[1] * 3)
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(alpha=0.2)
        _sample_note(ax, sample.kind, sample.doping, sample.T)
    return fig


def plot_ranking():
    """A physical crossing of lifetime curves, Eqs. (8), (35)."""
    fig, ax = _figure()
    dn = np.logspace(12, 17, 250)
    samples = [
        Sample(tau_n0=0.001, tau_p0=0.008, S_front=2, S_rear=2),
        Sample(tau_n0=0.002, tau_p0=0.0002, S_front=2, S_rear=2),
    ]
    for s, label in zip(
        samples, ("A: faster electron capture", "B: faster hole capture")
    ):
        ax.loglog(dn, component_lifetimes(dn, s)["effective"] * 1e6, label=label)
    ax.axvline(1e15, ls=":", color=".4", label="reporting point: 10$^{15}$ cm$^{-3}$")
    ax.set(
        xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
        ylabel="Effective lifetime (µs)",
        title="A ranking requires an injection level",
    )
    ax.legend(fontsize=8.5, loc="lower left", bbox_to_anchor=(0.02, 0.28))
    ax.grid(alpha=0.2)
    _sample_note(ax, s.kind, s.doping, s.T)
    return fig


def plot_voltage():
    """Implied-voltage curves and matched-generation operating points, Eqs. (43), (44)."""
    fig, axs = _figure(2)
    dn = np.logspace(12, 17, 220)
    for N in (1e15, 1e16):
        s = Sample(doping=N)
        axs[0].semilogx(
            dn,
            implied_voltage(dn, s) * 1e3,
            label=rf"$N_A={_scientific(N)}\ \mathrm{{cm}}^{{-3}}$",
        )
    rows = []
    for S in (1, 10, 100):
        s = Sample(S_front=S, S_rear=S)
        op = operating_point(s)
        rows.append((S, op))
        d = component_lifetimes(dn, s)
        G = dn / d["effective"]
        axs[1].semilogx(
            G, implied_voltage(dn, s) * 1e3, label=f"S = {S:g} cm/s per face"
        )
        axs[1].plot(
            op["G"], op["voltage"] * 1e3, "o", color=axs[1].lines[-1].get_color()
        )
    axs[0].set(
        xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
        ylabel="Implied voltage (mV)",
        title="Density maps to a quasi-Fermi splitting",
    )
    axs[1].axvline(rows[0][1]["G"], ls=":", color=".4")
    axs[1].set(
        xlabel="Absorbed generation G (cm$^{-3}$ s$^{-1}$)",
        ylabel="Implied voltage (mV)",
        title="Compare passivation at the same generation",
    )
    for ax in axs:
        ax.legend(fontsize=8.5, loc="upper left")
        ax.grid(alpha=0.2)
    _sample_note(
        axs[0], "p", r"$N_A$: values in legend", 300.0, xy=(0.96, 0.05), ha="right"
    )
    _sample_note(axs[1], s.kind, s.doping, s.T, xy=(0.96, 0.05), ha="right")
    return fig, rows


def plot_temperature_models():
    """Temperature changes ni and recombination coefficients, Eqs. (27), (19), (12)."""
    fig, axs = _figure(2)
    temp = np.linspace(275, 363, 120)
    axs[0].semilogy(temp, intrinsic_density(temp))
    axs[0].set(
        xlabel="Temperature (K)",
        ylabel="Intrinsic density (cm$^{-3}$)",
        title="Thermal carriers increase rapidly",
    )
    for dn in (1e14, 1e16):
        vals = [component_lifetimes(dn, Sample(T=float(T)))["intrinsic"] for T in temp]
        axs[1].plot(
            temp,
            np.asarray(vals) * 1e6,
            label=rf"$\Delta n={_scientific(dn)}\ \mathrm{{cm}}^{{-3}}$",
        )
    axs[1].set(
        xlabel="Temperature (K)",
        ylabel="Intrinsic lifetime (µs)",
        title="Intrinsic lifetime changes at fixed injection",
    )
    axs[1].legend(fontsize=9, loc="upper left")
    _sample_note(
        axs[1], "p", Sample().doping, "T: horizontal axis", xy=(0.96, 0.36), ha="right"
    )
    axs[0].text(
        0.04,
        0.95,
        "Intrinsic Si (undoped)\nT: horizontal axis",
        transform=axs[0].transAxes,
        va="top",
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.92, pad=3),
    )
    for ax in axs:
        ax.grid(alpha=0.2)
    return fig


def plot_bgn_sensitivity():
    """One-pass Schenk sensitivity, Eq. (31); not a coupled solver."""
    fig, axs = _figure(2)
    dn = np.logspace(12, 17, 180)
    s = Sample()
    n, p, n0, p0, ni = carrier_populations(dn, s.doping, s.kind, s.T)
    dEg = schenk_gap_narrowing(n, p, 0, s.doping, s.T)
    ni_eff = ni * np.exp(dEg / (2 * KB_EV * s.T))
    axs[0].semilogx(dn, dEg * 1e3)
    axs[0].set(
        ylabel="Gap reduction (meV)",
        title="Schenk evaluated at the unperturbed populations",
    )
    axs[1].semilogx(dn, implied_voltage(dn, s) * 1e3, label="unperturbed ni")
    axs[1].semilogx(
        dn,
        implied_voltage(dn, s, ni=ni_eff) * 1e3,
        "--",
        label="one-pass effective-ni illustration",
    )
    axs[1].set(
        ylabel="Implied voltage (mV)", title="The carrier-statistics convention matters"
    )
    axs[1].legend(fontsize=8)
    for ax in axs:
        ax.set_xlabel(r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)")
        ax.grid(alpha=0.2)
        _sample_note(ax, s.kind, s.doping, s.T, xy=(0.96, 0.05), ha="right")
    return fig


def plot_interface_physics():
    """Chemical and field-effect passivation, Eqs. (20), (22)."""
    fig, axs = _figure(2)
    E = np.linspace(-0.45, 0.45, 401)
    psi = np.linspace(-0.15, 0.15, 201)
    dn = 1e15
    for peak in (1e10, 1e11):
        Dit = interface_distribution(E, peak=peak)
        axs[0].plot(E, Dit, label=rf"${_scientific(peak)}$")
        S = interface_recombination(dn, psi, E, Dit) / dn
        axs[1].semilogy(psi, S, label=rf"${_scientific(peak)}$")
    axs[0].set(
        xlabel=r"Trap energy $E-E_i$ (eV)",
        ylabel=r"$D_{\mathrm{it}}$ ($\mathrm{cm}^{-2}\,\mathrm{eV}^{-1}$)",
        title="Chemical passivation removes states",
    )
    axs[1].set(
        xlabel="Prescribed surface potential (V)",
        ylabel="Effective surface velocity (cm/s)",
        title="Field effect changes which carriers reach them",
    )
    for ax in axs:
        ax.legend(
            fontsize=9,
            loc="upper left",
            title=r"Peak $D_{\mathrm{it}}$ ($\mathrm{cm}^{-2}\,\mathrm{eV}^{-1}$)",
            title_fontsize=9,
        )
        ax.grid(alpha=0.2)
    _sample_note(
        axs[1],
        "p",
        Sample().doping,
        300.0,
        xy=(0.96, 0.95),
        ha="right",
        va="top",
        extra=r"$\Delta n=10^{15}\ \mathrm{cm}^{-3}$",
    )
    return fig



def _temperature_colours(temperatures):
    """Map temperature in K from blue (low) to red (high) on a linear scale.

    Return a ScalarMappable shared by the curves and their colour bar.
    """
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    return ScalarMappable(
        norm=Normalize(min(temperatures), max(temperatures)), cmap="coolwarm"
    )


def plot_model_comparisons(data):
    """Doping, thickness and temperature sweeps using the same four channels.

    Each panel changes one input group, with all other inputs fixed. The
    optical and surface models are recalculated, rather than scaling τ_eff.
    Lifetime components combine through Eq. (37).
    """
    fig, axes = _figure(3, width=13, height=4.7, sharey=True)
    fig.supylabel("Effective lifetime (µs)")
    dn = data["delta_n_cm3"]
    sc = data["scenarios"]
    groups = [
        (
            ["doping_1e+15", "baseline", "doping_1e+16"],
            "N_dop",
            "Doping (n-type)",
            lambda v: rf"${_scientific(abs(v))}\ \mathrm{{cm}}^{{-3}}$",
        ),
        (
            ["thickness_100", "baseline", "thickness_300"],
            "W_um",
            "Thickness",
            lambda v: f"{v:g} µm",
        ),
        (
            ["temperature_275", "baseline", "temperature_350"],
            "T",
            "Temperature",
            lambda v: f"{v:g} K",
        ),
    ]
    for ax, (keys, parameter, title, fmt) in zip(axes, groups):
        temperatures = [sc[key]["parameters"]["T"] for key in keys]
        colours = _temperature_colours(temperatures) if parameter == "T" else None
        for key in keys:
            row = sc[key]
            style = (
                {"color": colours.to_rgba(row["parameters"]["T"])} if colours else {}
            )
            ax.loglog(
                dn,
                row["lifetimes_s"]["effective"] * 1e6,
                label=fmt(row["parameters"][parameter]),
                **style,
            )
        ax.set(xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)", title=title)
        if colours:
            fig.colorbar(
                colours,
                ax=ax,
                label="Temperature (K)",
                ticks=sorted(temperatures),
                pad=0.03,
                fraction=0.05,
            )
        else:
            ax.legend(fontsize=9, loc="upper right")
        ax.grid(alpha=0.2)
        baseline = sc["baseline"]["parameters"]
        density = (
            r"$N_D$: values in legend"
            if parameter == "N_dop"
            else abs(baseline["N_dop"])
        )
        temp = "T: colour bar" if parameter == "T" else baseline["T"]
        _sample_note(
            ax, "n" if baseline["N_dop"] < 0 else "p", density, temp, fontsize=8.5
        )
    return fig


def plot_injection_channels(data):
    """Injection trends of all four channels for the baseline model sample.

    Eq. (39) relates slope to recombination order. Dashed curves
    expose simple constant-coefficient limits, not fits to the curves.
    """
    fig, axes = _figure(2, width=10.5, nrows=2, height=7.2)
    dn = data["delta_n_cm3"]
    baseline = data["scenarios"]["baseline"]
    d = baseline["lifetimes_s"]
    params = baseline["parameters"]
    for ax, key, title in zip(
        axes.flat,
        ["radiative", "auger", "srh", "surface"],
        [
            "Radiative: pair encounters",
            "Auger: a third carrier",
            "Bulk SRH: sequential capture",
            "Surface: fixed S_eff (flat by construction)",
        ],
    ):
        ax.loglog(dn, d[key] * 1e6, lw=2, label="Lifetime model calculation")
        if key in ("radiative", "auger"):
            power = -1 if key == "radiative" else -2
            ids = dn >= 1e16
            ax.loglog(
                dn[ids],
                d[key][-1] * 1e6 * (dn[ids] / dn[-1]) ** power,
                "--",
                color=".4",
                label=f"High-injection slope {power}",
            )
        ax.axvline(
            abs(params["N_dop"]),
            color=".5",
            ls=":",
            label=r"$\Delta n=|N_{\mathrm{dop}}|$",
        )
        ax.set(
            xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
            ylabel="Lifetime (µs)",
            title=title,
        )
        if key == "radiative":
            ax.set_ylim(top=ax.get_ylim()[1] * 1.8)
        ax.legend(fontsize=8, loc="upper right" if key != "srh" else "upper left")
        ax.grid(alpha=0.2)
        _sample_note(
            ax,
            "n" if params["N_dop"] < 0 else "p",
            abs(params["N_dop"]),
            params["T"],
            xy=(0.96, 0.05) if key == "srh" else (0.04, 0.05),
            ha="right" if key == "srh" else "left",
        )
    return fig


def plot_recombination_coefficients():
    """Show temperature and carrier screening in Eqs. (12), (19).

    Coefficients are evaluated for a symmetric plasma n=p; this isolates the
    coefficient model from equilibrium statistics and sample electrostatics.
    """
    from matplotlib.lines import Line2D

    fig, axes = _figure(2)
    density = np.logspace(13, 19, 200)
    temperatures = (275.0, 300.0, 350.0)
    colours = _temperature_colours(temperatures)
    for T in temperatures:
        cn, cp = auger_coefficients(density / 2, density / 2, T)
        colour = colours.to_rgba(T)
        axes[0].loglog(density, cn, color=colour, lw=2)
        axes[0].loglog(density, cp, "--", color=colour, lw=2)
        B = radiative_coefficient(T) * b_relative(density / 2, density / 2, T)
        axes[1].loglog(density, B, color=colour, lw=2)
    axes[0].set(
        ylabel="Auger coefficient (cm$^6$/s)",
        title="Black: screened and enhanced channels",
    )
    axes[1].set(
        ylabel="B before optical recycling (cm$^3$/s)",
        title="Nguyen × Altermatt density correction",
    )
    axes[0].legend(
        handles=[
            Line2D([], [], color=".25", label="$C_n$"),
            Line2D([], [], color=".25", ls="--", label="$C_p$"),
        ],
        fontsize=9,
        ncol=2,
    )
    fig.colorbar(
        colours,
        ax=list(axes),
        label="Temperature (K)",
        ticks=temperatures,
        pad=0.03,
        fraction=0.025,
    )
    for ax in axes:
        ax.set_xlabel("Total carrier density n + p (cm$^{-3}$)")
        ax.grid(alpha=0.2)
    axes[1].text(
        0.04,
        0.06,
        "Symmetric plasma: n = p\nDoping / type not specified\nT: colour bar",
        transform=axes[1].transAxes,
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.92, pad=3),
    )
    return fig


def print_voltage_examples(examples):
    """Print the operating points from plot_voltage, Eqs. (43), (44).

    examples : iterable of (S, point) pairs; S is the velocity per face in
        cm/s, point['delta_n'] is in cm^-3 and point['voltage'] is in V.
    Returns None. Only the displayed voltage is converted from V to mV.
    """
    for velocity, point in examples:
        print(
            f"S={velocity:3g} cm/s: delta_n={point['delta_n']:.3e} cm^-3, "
            f"implied voltage={point['voltage']*1e3:.1f} mV"
        )


def plot_model_components(data):
    """The four lifetime contributions and their harmonic sum, Eq. (37).

    data comes from model_examples(). These are forward-model curves, not
    measurements.
    """
    fig, ax = _figure(width=8.0)
    fig.set_size_inches(8.0, 5.8)
    d = data["scenarios"]["baseline"]["lifetimes_s"]
    dn = data["delta_n_cm3"]
    for key, label, color in [
        ("srh", "Bulk SRH", "C0"),
        ("radiative", "Radiative", "C1"),
        ("auger", "Auger (Black model)", "C2"),
        ("surface", "Surface", "C3"),
        ("effective", "Effective", "k"),
    ]:
        ax.loglog(
            dn,
            d[key] * 1e6,
            label=label,
            color=color,
            lw=2.8 if key == "effective" else 1.8,
        )
    ax.loglog(
        dn,
        d["intrinsic"] * 1e6,
        "--",
        color=".5",
        lw=1.2,
        label="Intrinsic (radiative + Auger)",
    )
    ax.set(
        xlabel=r"Excess density $\Delta n$ ($\mathrm{cm}^{-3}$)",
        ylabel="Lifetime (µs)",
        title="Four recombination channels and effective lifetime",
    )
    ax.legend(fontsize=9, ncol=2, loc="lower left")
    ax.grid(alpha=0.2)
    params = data["scenarios"]["baseline"]["parameters"]
    _sample_note(
        ax,
        "n" if params["N_dop"] < 0 else "p",
        abs(params["N_dop"]),
        params["T"],
        xy=(0.04, 0.19),
        extra=f"W = {params['W_um']:g} µm",
    )
    return fig


# ---------------------------------------------------------------------------
# Model curves for the combined-lifetime sections (6 and 7)
# ---------------------------------------------------------------------------
def model_examples(delta_n=None):
    """Forward-model lifetime curves for one baseline sample and its sweeps.

    Returns {"delta_n_cm3": array, "scenarios": {label: {"parameters": {...},
    "lifetimes_s": {srh, radiative, auger, surface, intrinsic, effective}}}}.
    The baseline is a passivated n-type wafer (3e15 cm^-3, 150 um, 300 K,
    S = 10 cm/s per face, generic deep level). Doping, thickness and
    temperature are then changed one at a time. Every curve is computed
    from component_lifetimes, so the rates combine through Eq. (37).
    These are model curves, not measurements. The surface term uses a
    fixed S_eff; its injection dependence is the subject of Sec. 4.4.
    """
    dn = np.logspace(13, 17, 65) if delta_n is None else np.asarray(delta_n, float)
    base = Sample(doping=3e15, kind="n", T=300.0, W_um=150.0, S_front=10.0, S_rear=10.0)
    changes = [("baseline", {})]
    changes += [(f"doping_{N:.0e}", dict(doping=N)) for N in (1e15, 1e16)]
    changes += [(f"thickness_{W:g}", dict(W_um=W)) for W in (100.0, 300.0)]
    changes += [(f"temperature_{T:g}", dict(T=T)) for T in (275.0, 350.0)]
    scenarios = {}
    for label, change in changes:
        s = replace(base, **change)
        parts = component_lifetimes(dn, s)
        scenarios[label] = dict(
            parameters=dict(
                N_dop=-s.doping if s.kind == "n" else s.doping,
                T=s.T, W_um=s.W_um, S_front=s.S_front, S_rear=s.S_rear),
            lifetimes_s={k: np.asarray(v) for k, v in parts.items() if k != "bulk"},
        )
    return dict(delta_n_cm3=dn, scenarios=scenarios)


def plot_recombination_mechanisms():
    """Band-diagram schematic of the four recombination mechanisms.

    Radiative, Auger (eeh and ehh), bulk SRH and surface SRH. Schematic only:
    the radiative panel omits the phonon that conserves crystal momentum in
    silicon. Arrows are carrier transitions; the wavy line is a photon.
    """
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.9), constrained_layout=True)
    ec, ev = 1.0, 0.0
    blue, red, grey = "#1f6fb4", "#c0392b", "#555555"

    def bands(ax, title, x0=0.0, x1=1.0, surface=False):
        ax.set(xlim=(x0, x1), ylim=(-0.45, 1.75))
        ax.axis("off")
        if surface:
            xs = np.linspace(x0, x1, 80)
            bend = 0.28 * np.exp(-((xs - x0) / 0.22) ** 2)
            ax.plot(xs, ec + bend, color=grey, lw=2)
            ax.plot(xs, ev + bend, color=grey, lw=2)
        else:
            ax.plot([x0, x1], [ec, ec], color=grey, lw=2)
            ax.plot([x0, x1], [ev, ev], color=grey, lw=2)
        ax.text(x1, ec + 0.07, "$E_c$", ha="right", fontsize=10)
        ax.text(x1, ev - 0.17, "$E_v$", ha="right", fontsize=10)
        ax.set_title(title, fontsize=11)

    def arrow(ax, x, y0, y1, color, style="-|>", ls="-"):
        ax.annotate("", xy=(x, y1), xytext=(x, y0),
                    arrowprops=dict(arrowstyle=style, color=color, lw=2, ls=ls))

    def electron(ax, x, y, filled=True):
        ax.plot([x], [y], "o", ms=9, mfc=blue if filled else "white", mec=blue)

    def hole(ax, x, y):
        ax.plot([x], [y], "o", ms=9, mfc="white", mec=red, mew=2)

    def photon(ax, x, y0, y1):
        ys = np.linspace(y0, y1, 60)
        ax.plot(x + 0.035 * np.sin(ys * 55), ys, color="#d4a017", lw=2)

    # (a) radiative
    ax = axes[0]
    bands(ax, "Radiative")
    electron(ax, 0.40, ec + 0.12)
    hole(ax, 0.40, ev - 0.12)
    arrow(ax, 0.40, ec + 0.05, ev + 0.05, blue)
    photon(ax, 0.68, 0.85, 0.15)
    ax.text(0.50, 0.50, "photon", fontsize=9, color="#9a7410", ha="left")
    # (b) Auger
    ax = axes[1]
    bands(ax, "Auger (eeh)")
    electron(ax, 0.30, ec + 0.12)
    electron(ax, 0.55, ec + 0.12)
    hole(ax, 0.30, ev - 0.12)
    arrow(ax, 0.30, ec + 0.05, ev + 0.05, blue)
    arrow(ax, 0.55, ec + 0.10, ec + 0.55, blue)
    arrow(ax, 0.70, ec + 0.55, ec + 0.12, blue, ls=":")
    ax.text(0.40, ec + 0.62, "energy to a third carrier, then heat", fontsize=8.5)
    # (c) bulk SRH
    ax = axes[2]
    bands(ax, "Bulk SRH")
    ax.plot([0.25, 0.75], [0.5, 0.5], color="#7d3c98", lw=2.5, ls="--")
    ax.text(0.77, 0.5, "$E_t$", va="center", fontsize=10, color="#7d3c98")
    electron(ax, 0.38, ec + 0.12)
    hole(ax, 0.62, ev - 0.12)
    arrow(ax, 0.38, ec + 0.05, 0.55, blue)
    arrow(ax, 0.62, ev - 0.05, 0.45, red)
    ax.text(0.5, -0.36, "electron captured, then hole captured", ha="center", fontsize=8.5)
    # (d) surface SRH
    ax = axes[3]
    bands(ax, "Surface SRH", surface=True)
    for y in np.linspace(0.38, 0.62, 5):
        ax.plot([0.06], [y], "s", ms=4, color="#7d3c98")
    ax.text(0.04, -0.36, "interface states (squares)", fontsize=8.5, color="#7d3c98")
    electron(ax, 0.55, ec + 0.12)
    hole(ax, 0.55, ev - 0.12)
    ax.annotate("", xy=(0.12, 0.62), xytext=(0.55, 1.05),
                arrowprops=dict(arrowstyle="-|>", color=blue, lw=2, connectionstyle="arc3,rad=0.15"))
    ax.annotate("", xy=(0.12, 0.38), xytext=(0.55, -0.05),
                arrowprops=dict(arrowstyle="-|>", color=red, lw=2, connectionstyle="arc3,rad=-0.15"))
    ax.text(0.04, 1.55, "surface", fontsize=8.5, color=grey)
    return fig




# ===========================================================================
# PART II - THE MEASUREMENT (Eqs. 45-57)
# ===========================================================================
# Dependencies: carrier_populations, intrinsic_density and the validators
# above are shared with Part I. Lengths are cm unless a name or docstring
# says mm.


def _positive(value, name, zero=False):
    """Validate finite positive (or nonnegative) scalar/array inputs."""
    a = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(a < 0 if zero else a <= 0):
        raise ValueError(f'{name} must be finite and {"nonnegative" if zero else "positive"}')
    return a


def equilibrium_carriers(doping_cm3, sample_type='n', ni_cm3=9.65e9):
    """Return equilibrium n0, p0 (cm^-3) using an explicit ni.

    Densities in cm^-3. The explicit ni is retained for callers comparing
    equilibrium statistics. Measurement calculations call carrier_populations
    with their actual T, using the temperature-dependent ni instead.
    """
    _, _, n0, p0, _ = carrier_populations(0., doping_cm3, sample_type, ni=ni_cm3)
    return n0, p0


def sheet_conductance(dn_cm3, doping_cm3, W_cm, sample_type='n', T_K=300.):
    """Total sheet conductance in S, Eq. (49), including dark carriers.

    Mobilities depend on total n,p, doping and T. Conductivity is assumed
    uniform through W; conducting coatings are outside this bulk model.
    """
    dn = _positive(dn_cm3, 'excess density', zero=True)
    W = float(_positive(W_cm, 'thickness'))
    n,p,n0,p0,_ = carrier_populations(dn, doping_cm3, sample_type, T_K)
    nd,na = (doping_cm3,0.) if sample_type == 'n' else (0.,doping_cm3)
    mn,mp = klaassen_unified_mobility_si(T_K,nd,na,n,p)
    return Q*W*(mn*n+mp*p)


def doping_from_resistivity(rho_ohm_cm, sample_type='n', T_K=300.):
    """Infer uncompensated doping (cm^-3) from dark rho using Eq. (49).

    Uses log-space bracketing over 1e10--1e19 cm^-3; fails if unbracketed.
    A conducting surface or compensation invalidates this inference.
    """
    rho = float(_positive(rho_ohm_cm, 'resistivity'))
    root = brentq(lambda x: 1/sheet_conductance(0.,10**x,1.,sample_type,T_K)-rho,10,19)
    return 10**root


def excess_density(total_sheet_S, doping_cm3, W_cm, sample_type='n', T_K=300.):
    """Invert Eq. (49) by bracketed root finding; return dn in cm^-3.

    Nonfinite or sub-dark conductances return NaN, not a clipped positive
    density. Exact dark returns zero. Roots above 1e19 are not extrapolated.
    """
    values = np.asarray(total_sheet_S,dtype=float)
    dark = float(sheet_conductance(0.,doping_cm3,W_cm,sample_type,T_K))
    out = np.full(values.shape,np.nan)
    for index in np.ndindex(values.shape):
        target=values[index]
        if not np.isfinite(target) or target < dark:
            continue
        if target == dark:
            out[index]=0.; continue
        def residual(x):
            """Conductance residual at log10(1+dn)."""
            return sheet_conductance(10**x-1,doping_cm3,W_cm,sample_type,T_K)-target
        out[index]=10**brentq(residual,0,19,xtol=1e-11)-1
    return out


def layer_sensitivity(z_mm, lambda_mm):
    """Dimensionless S'(z), Eq. (46), Black & Macdonald 2019 Eq. (4)."""
    z = _positive(z_mm,'depth',zero=True)
    lam = float(_positive(lambda_mm,'lambda'))
    return (1+z/lam)**-2/(1+(z/(2*lam))**2)


def average_sensitivity(W_mm, lift_mm=0., lambda_mm=3.9):
    """Uniform-slab sensitivity, Eq. (47); all geometry inputs in mm.

    Zero lift is the fitted stage reference, not the physical wire surface.
    This is the absolute thin-sheet convention, S'(0)=1.
    """
    W=float(_positive(W_mm,'W')); z=float(_positive(lift_mm,'lift',zero=True))
    return quad(lambda u: float(layer_sensitivity(u,lambda_mm)),z,z+W)[0]/W


def voltage_to_sheet(V_sample, V_air, coefficients, W_mm, lift_mm=0., lambda_mm=3.9):
    """Voltage to actual total sheet S using Eqs. (45)--(47).

    coefficients=(A_c,B_c,C_c), in S/V^2,S/V,S; x=V_sample-V_air.
    Must use coefficients fitted in the same absolute sensitivity convention.
    Do not apply this correction a second time to an already-corrected export.
    """
    x=np.asarray(V_sample,dtype=float)-V_air
    return np.polyval(coefficients,x)/average_sensitivity(W_mm,lift_mm,lambda_mm)


def fit_coil_calibration(rows=None, air=None, resistance=None, W_mm=.525, max_training_layers=2):
    """Fit lambda and A_c/B_c/C_c, Eq. (48), from known wafer stacks.

    rows: (date, tuple of wafer IDs bottom-to-top, voltage replicates).
    air: date->V_air; resistance: wafer ID->sheet ohms. Defaults are the
    64 measured configurations embedded below (CALIBRATION_ROWS). Singles
    and pairs fit by default; triples/quadruples are held-out diagnostics.
    For each lambda, re-fit a quadratic by relative least squares. Prefix
    inversion reconstructs each layer's conductance; it assumes unchanged
    lower layers, no extra gaps, and stable calibration between dates.
    Returns parameters, profile, per-configuration table and grouped errors.
    """
    rows = CALIBRATION_ROWS if rows is None else rows
    air = CALIBRATION_AIR if air is None else air
    resistance = CALIBRATION_RESISTANCE if resistance is None else resistance
    stacks=[tuple(r[1]) for r in rows]
    count=np.array([len(s) for s in stacks])
    x=np.array([np.mean(r[2])-air[r[0]] for r in rows])
    X=np.column_stack([x*x,x,np.ones_like(x)])
    H=np.zeros((len(rows),int(count.max())))
    for k,s in enumerate(stacks):
        for j,key in enumerate(s): H[k,j]=1/resistance[key]
    train=count<=max_training_layers
    if train.sum()<4: raise ValueError('At least four calibration configurations required')
    def profile(lam):
        """Relative residual objective after refitting polynomial at lambda."""
        sens=np.array([average_sensitivity(W_mm,i*W_mm,lam) for i in range(H.shape[1])])
        target=H@sens
        abc=np.linalg.lstsq(X[train]/target[train,None],np.ones(train.sum()),rcond=None)[0]
        return np.mean((X[train]@abc/target[train]-1)**2),abc,target,sens
    grid=np.geomspace(.05,100.,180)
    objective=np.array([profile(l)[0] for l in grid])
    k=int(np.argmin(objective))
    if k in (0,len(grid)-1): raise ValueError('Coil attenuation optimum is at search boundary')
    opt=minimize_scalar(lambda l:profile(l)[0],bounds=(grid[k-1],grid[k+1]),method='bounded')
    loss,abc,target,sens=profile(opt.x)
    apparent=X@abc; lookup=dict(zip(stacks,apparent)); corrected=[]
    for s in stacks:
        try: increments=np.diff([0.]+[lookup[s[:i]] for i in range(1,len(s)+1)])
        except KeyError as exc: raise ValueError('Missing measured stack prefix') from exc
        corrected.append(np.sum(increments/sens[:len(s)]))
    standard=H.sum(axis=1); corrected=np.array(corrected)
    table=pd.DataFrame({'layers':count,'training':train,'standard_S':standard,
                        'apparent_S':apparent,'corrected_S':corrected,
                        'error_percent':100*(corrected/standard-1)})
    summary=table.groupby('layers').agg(N=('layers','size'),
        mean_error_percent=('error_percent','mean'),
        mean_absolute_error_percent=('error_percent',lambda v:np.mean(abs(v))))
    return dict(lambda_mm=float(opt.x),coefficients=abc,loss=loss,table=table,
                summary=summary,grid_mm=grid,objective=objective,voltage=x,target_S=target)


def generation_rate(V_ref, W_cm, C_opt, ref_V_per_sun, J_ref_A_cm2=.038, V_offset=0.):
    """G in cm^-3 s^-1, Eq. (51), from relative reference illumination.

    C_opt is an effective scale, not strictly an absorptance (can exceed 1).
    ref_V_per_sun is specific to reference aperture/gain; J_ref=.038 A/cm^2
    is the Sinton generation-current convention. No negative clipping.
    """
    for value,name in [(W_cm,'W'),(C_opt,'C_opt'),(ref_V_per_sun,'reference calibration'),(J_ref_A_cm2,'J_ref')]:
        _positive(value,name)
    return (np.asarray(V_ref)-V_offset)/ref_V_per_sun*C_opt*J_ref_A_cm2/(Q*W_cm)


def reference_calibration(V_ref, known_suns):
    """Fit V_ref=offset+gain*suns to independent known light levels.

    Returns (gain V/sun, offset V, RMS V). At least three finite distinct
    levels are required. A fitted gain is not a traceable calibration unless
    the independent light standard, spectrum and aperture are controlled.
    """
    x=np.asarray(known_suns,float); y=np.asarray(V_ref,float)
    if x.shape!=y.shape or x.ndim!=1 or len(x)<3 or len(np.unique(x))<3 or not np.isfinite(x+y).all():
        raise ValueError('Need three or more paired finite distinct light levels')
    gain,offset=np.polyfit(x,y,1)
    _positive(gain,'reference gain')
    return gain,offset,float(np.sqrt(np.mean((y-gain*x-offset)**2)))


def carrier_derivative(time_s, dn_cm3, window=5, polyorder=2):
    """Return (dn/dt, valid_mask), using local polynomial ln(dn) fits.

    Savitzky--Golay on uniform time only. Fit each contiguous positive,
    finite run independently; mask half-window edges. Nonuniform time,
    duplicate timestamps and invalid window choices raise ValueError.
    No log-clipping, interpolation across gaps, or automatic window switch.
    """
    t=np.asarray(time_s,float); n=np.asarray(dn_cm3,float)
    if t.ndim!=1 or t.shape!=n.shape or len(t)<3 or not np.isfinite(t).all():
        raise ValueError('Need matching 1-D time and density arrays')
    dt=np.diff(t)
    if np.any(dt<=0) or not np.allclose(dt,np.median(dt),rtol=1e-4,atol=1e-12):
        raise ValueError('Savitzky-Golay requires strictly increasing uniformly spaced times')
    if window%2!=1 or window<3 or polyorder<1 or polyorder>=window:
        raise ValueError('Use odd window >=3 and 1 <= polyorder < window')
    good=np.isfinite(n)&(n>0); derivative=np.full_like(n,np.nan)
    boundaries=np.flatnonzero(np.diff(np.r_[False,good,False]))
    half=window//2
    for start,end in zip(boundaries[::2],boundaries[1::2]):
        if end-start<window: continue
        d=savgol_filter(np.log(n[start:end]),window,polyorder,deriv=1,delta=np.median(dt))
        derivative[start+half:end-half]=n[start+half:end-half]*d[half:-half]
    return derivative,np.isfinite(derivative)


def lifetime_from_signals(time_s, dn_cm3, G_cm3_s, mode='generalized', window=5):
    """Return lifetime/diagnostics from Eqs. (52)--(54).

    valid additionally restricts to the post-peak, decreasing-density branch.
    qss_error=|dn/dt|/G is exact relative QSS error vs generalized balance;
    transient_error=G/(-dn/dt) is exact relative transient error on decay.
    Positive lifetime alone does not imply physical validity.
    """
    n=np.asarray(dn_cm3,float); g=np.broadcast_to(np.asarray(G_cm3_s,float),n.shape)
    derivative,valid_d=carrier_derivative(time_s,n,window)
    if mode not in ('generalized','qss','transient'): raise ValueError('Unknown analysis mode')
    rate={'generalized':g-derivative,'qss':g,'transient':-derivative}[mode]
    valid=valid_d & np.isfinite(g)&(g>=0)&(rate>0)&(derivative<0)
    if np.isfinite(n).any(): valid &= np.arange(n.size)>np.nanargmax(n)
    tau=np.full_like(n,np.nan)
    np.divide(n,rate,out=tau,where=valid)
    with np.errstate(divide='ignore',invalid='ignore'):
        qss=np.abs(derivative)/g; transient=g/(-derivative)
    return dict(dn=n,G=g,derivative=derivative,tau=tau,valid=valid,
                qss_error=qss,transient_error=transient,time=np.asarray(time_s),mode=mode)


def load_sinton(path=None, *, text=None):
    """Parse a tab-delimited Sinton export; returns metadata, DataFrame.

    Header is found by Time_(s); metadata are read only before that header.
    Numerical NaN remains NaN. Data are not sorted or silently dropped.
    """
    if (path is None)==(text is None): raise ValueError('Provide exactly one of path or text')
    text=Path(path).read_text(encoding='utf-8-sig') if text is None else text
    lines=text.splitlines()
    try: header=next(i for i,line in enumerate(lines) if line.startswith('Time_(s)\t'))
    except StopIteration as exc: raise ValueError('Missing Time_(s) header') from exc
    meta={}
    for line in lines[:header]:
        parts=line.split('\t',1)
        if len(parts)==2: meta[parts[0]]=parts[1].strip()
    frame=pd.read_csv(io.StringIO('\n'.join(lines[header:])),sep='\t')
    for col in frame: frame[col]=pd.to_numeric(frame[col],errors='raise')
    if np.any(np.diff(frame['Time_(s)'])<=0): raise ValueError('Time must strictly increase')
    return meta,frame


def example_measurement(name='3.txt'):
    """Load an exact embedded example export, verified by SHA256.

    Base85/zlib stores data only, never executable source; this keeps the
    requested two-file package portable. See SOURCE_SHA256 for provenance.
    """
    raw=zlib.decompress(base64.b85decode(EMBEDDED_EXPORTS[name]))
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA256[name]: raise ValueError('Data checksum mismatch')
    return load_sinton(text=raw.decode('utf-8-sig'))


def analyse_measurement(measurement, C_opt=None, ref_V_per_sun=.0113756, window=5, mode=None):
    """Reduce one export using its already-calibrated excess sheet column.

    Adds 1/R_sheet,dark to Conductivity_(Siemens), (the export column is named Conductivity but holds excess sheet
    conductance).
    Infers doping consistently from measured rho. Does NOT reapply a coil
    calibration to exports without matching records. The default reference gain
    is a value carried over from earlier analysis, not independently verified
    instrument metadata.
    """
    meta,df=measurement
    st={'n-type':'n','p-type':'p'}[meta['Sample_Type'].lower()]
    T=float(meta['Temperature_(C)'])+273.15; W=float(meta['Thickness_(cm)'])
    rho=float(meta['Measured_Resistivity_(ohm-cm)'])
    N=doping_from_resistivity(rho,st,T)
    sigma=df['Conductivity_(Siemens)'].to_numpy()+1/float(meta['Total_Sheet_Resistance_(ohm/sq)'])
    n=excess_density(sigma,N,W,st,T)
    c=float(meta['Optical_Constant']) if C_opt is None else C_opt
    g=generation_rate(df['Reference_Voltage_(V)'].to_numpy(),W,c,ref_V_per_sun)
    mode=mode or ('transient' if meta['Analysis_Mode'].lower()=='transient' else 'generalized')
    result=lifetime_from_signals(df['Time_(s)'].to_numpy(),n,g,mode,window)
    result.update(meta=meta,frame=df,doping=N,T_K=T,W_cm=W,C_opt=c,ref_V_per_sun=ref_V_per_sun)
    return result


def _curve(n,y,mask):
    """Sorted unique positive curve; no extrapolation or duplicate averaging."""
    good=np.asarray(mask)&np.isfinite(n)&np.isfinite(y)&(n>0)&(y>0)
    x=np.asarray(n)[good]; v=np.asarray(y)[good]; order=np.argsort(x)
    x=x[order];v=v[order]
    if len(x)<5 or np.any(np.diff(x)<=0): raise ValueError('Need >=5 unique valid injection points')
    return x,v


def fit_optical_constant(transient, generalized, injection_range=(2e15,7e15),
                         bounds=(.3,3.), max_transient_error=.05):
    """Fit C_opt on fixed trustworthy overlap, Eq. (55), in log lifetime.

    Require residual generation <=5% of decay in transient data, using its
    stated optical scale; require positive generalized recombination at BOTH
    C_opt bounds. Compare on a fixed log-injection grid, never extrapolate.
    The returned RMS is a mismatch diagnostic, not a statistical uncertainty.
    This cannot separate reference gain error from an optical scale error.
    """
    t,g=transient,generalized
    if not 0<bounds[0]<bounds[1]: raise ValueError('Invalid positive C_opt bounds')
    mt=t['valid']&(t['transient_error']<=max_transient_error)
    xt,yt=_curve(t['dn'],t['tau'],mt)
    Gunit=g['G']/g['C_opt']
    mg=g['valid']&(bounds[0]*Gunit-g['derivative']>0)&(bounds[1]*Gunit-g['derivative']>0)
    xg,_=_curve(g['dn'],g['tau'],mg)
    lo=max(injection_range[0],xt.min(),xg.min()); hi=min(injection_range[1],xt.max(),xg.max())
    if lo>=hi: raise ValueError('No trustworthy common injection range')
    if np.sum((xt>=lo)&(xt<=hi))<5 or np.sum((xg>=lo)&(xg<=hi))<5:
        raise ValueError('Fewer than five measured points per trace in fitting window')
    grid=np.geomspace(lo,hi,100); target=np.interp(np.log10(grid),np.log10(xt),np.log10(yt))
    def objective(c):
        """Mean squared log10 lifetime mismatch on a fixed grid."""
        tau=np.full_like(g['dn'],np.nan)
        np.divide(g['dn'],c*Gunit-g['derivative'],out=tau,where=mg)
        xx,yy=_curve(g['dn'],tau,mg)
        residual=np.interp(np.log10(grid),np.log10(xx),np.log10(yy))-target
        return np.mean(residual**2)
    fit=minimize_scalar(objective,bounds=bounds,method='bounded')
    if not fit.success or min(fit.x-bounds[0],bounds[1]-fit.x)<.001*(bounds[1]-bounds[0]):
        raise ValueError('Optical fit failed or reached a bound')
    scan=np.linspace(*bounds,100)
    return dict(C_opt=float(fit.x),rms_decades=float(np.sqrt(fit.fun)),
                injection_range=(lo,hi),grid=grid,scan=scan,
                objective=np.array([objective(c) for c in scan]))


def lifetime_at_injection(result, injection_cm3, extra_mask=None):
    """Log-log interpolate a valid curve at dn; reject extrapolation."""
    mask=result['valid'].copy()
    if extra_mask is not None: mask &= extra_mask
    x,y=_curve(result['dn'],result['tau'],mask)
    if not x.min()<=injection_cm3<=x.max(): raise ValueError('Requested injection outside valid data')
    return 10**np.interp(np.log10(injection_cm3),np.log10(x),np.log10(y))


def plot_instrument():
    """Original schematic of measured channels and inferred quantities."""
    fig,ax=plt.subplots(figsize=(11,4.6),constrained_layout=True)
    ax.set(xlim=(0,11),ylim=(0,5));ax.axis('off')
    boxes=[(1.1,4.4,'Flash lamp + filter'),(1.1,2.7,'Wafer on temperature stage'),
           (1.1,.9,'RF coil / electronics'),(6.,4.4,'Reference cell'),
           (6.,2.7,'Reference voltage → G(t)'),(6.,.9,'Coil voltage → sheet S → Δn(t)'),
           (9.1,1.8,'Carrier balance\n→ effective lifetime')]
    for x,y,label in boxes:
        ax.text(x,y,label,ha='center',va='center',fontsize=10,
                bbox=dict(boxstyle='round,pad=.5',fc='#eef4f8',ec='#42677e'))
    for start,end in [((1.1,4.05),(1.1,3.1)),((2.8,4.4),(4.9,4.4)),((6,4.05),(6,3.1)),
                      ((1.1,2.35),(1.1,1.3)),((2.5,.9),(4,.9)),((7.4,2.7),(8.6,2.1)),((7.6,.9),(8.6,1.5))]:
        ax.annotate('',xy=end,xytext=start,arrowprops=dict(arrowstyle='->',color='#42677e',lw=1.5))
    ax.set_title('What is measured, and what is inferred',loc='left')
    return fig


def plot_sensitivity(lambda_mm=3.9):
    """Thickness and lift-off effects, Eqs. (46)--(47); model curves."""
    fig,axes=plt.subplots(1,2,figsize=(10,3.6),constrained_layout=True)
    z=np.linspace(0,3,150);axes[0].plot(z,layer_sensitivity(z,lambda_mm))
    axes[0].set(xlabel='Layer distance from stage reference (mm)',ylabel="Relative sensitivity S′",title='Thin conducting layer')
    w=np.linspace(.05,2.1,100)
    for lift in [0,.2,.5]: axes[1].plot(w,[average_sensitivity(v,lift,lambda_mm) for v in w],label=f'Lift-off {lift:g} mm')
    axes[1].set(xlabel='Wafer thickness (mm)',ylabel='Mean sensitivity',title=f'Uniform slab; λ = {lambda_mm:.2f} mm')
    axes[1].legend()
    for ax in axes: ax.grid(alpha=.2)
    return fig


def plot_coil_calibration(result):
    """Measured stack calibration and re-fitted objective profile."""
    fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
    d=result['table'];x=result['voltage'];grid=np.linspace(0,x.max(),100)
    axes[0].scatter(x,result['target_S'],s=15,c=d['layers'],cmap='viridis')
    axes[0].plot(grid,np.polyval(result['coefficients'],grid),'k--')
    axes[0].set(xlabel='Sample minus air voltage (V)',ylabel='Sensitivity-weighted sheet S (S)',title='Known standards and quadratic fit')
    axes[1].loglog(result['grid_mm'],result['objective'])
    axes[1].axvline(result['lambda_mm'],ls='--',color='k')
    axes[1].set(xlabel='λ (mm)',ylabel='Mean squared relative residual',title='Polynomial re-fitted at each λ')
    for training,label,marker in [(True,'Fit: singles / pairs','o'),(False,'Held out: triples / quads','^')]:
        sub=d[d.training==training];axes[2].scatter(sub.standard_S,sub.error_percent,label=label,marker=marker,s=18)
    axes[2].axhline(0,color='k',lw=.8);axes[2].legend(fontsize=8)
    axes[2].set(xlabel='Known total sheet S (S)',ylabel='Recovered conductance error (%)',title='Measured-prefix reconstruction')
    return fig


def plot_mode_limits():
    """Exact constant-lifetime exponential-flash example; synthetic, Eq. 53.

    Starts at dn=0; d(dn)/dt=G0 exp(-t/tf)-dn/tau is solved analytically.
    Illustrates mode error after peak without derivative-estimation noise.
    """
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
    tf=.002;t=np.linspace(1e-6,.018,1000);g0=1e19;g=g0*np.exp(-t/tf)
    for tau in [50e-6,200e-6,1e-3]:
        n=g0*(np.exp(-t/tf)-np.exp(-t/tau))/(1/tau-1/tf)
        dn=g-n/tau;decay=dn<0
        axes[0].plot(t*1e3,n/1e15,label=f'τ = {tau*1e6:g} µs')
        axes[1].plot(t[decay]*1e3,(n/g)[decay]/tau)
        axes[2].plot(t[decay]*1e3,(-n/dn)[decay]/tau)
    axes[0].legend(fontsize=9);axes[0].set(ylabel='Excess density (10¹⁵ cm⁻³)',title='Carrier storage under a 2 ms flash')
    for ax,label in zip(axes[1:],['QSS estimate / true lifetime','Transient estimate / true lifetime']):
        ax.axhline(1,color='k',ls='--');ax.set(ylim=(0,4),ylabel=label)
    for ax in axes: ax.set_xlabel('Time (ms)');ax.grid(alpha=.2)
    return fig


def plot_raw_measurements(transient,generalized):
    """Raw channels and mobility-reconstructed density from two exports."""
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
    for r,label in [(transient,'3.txt: short flash'),(generalized,'4.txt: full flash')]:
        df=r['frame'];t=r['time']*1e3
        axes[0].plot(t,df['Photoconductance_Signal_(V)'],label=label)
        axes[1].semilogy(t,np.where(df['Reference_Voltage_(V)']>0,df['Reference_Voltage_(V)'],np.nan))
        axes[2].semilogy(t,r['dn'])
    for ax in axes: ax.set_xlabel('Time (ms)');ax.grid(alpha=.2)
    axes[0].set_ylabel('Coil channel (V)');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Reference channel (V)');axes[2].set_ylabel('Inferred Δn (cm⁻³)')
    return fig


def plot_optical_fit(t,g,fit):
    """Measured mode matching and diagnostic objective; not a precision claim."""
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
    c=fit['C_opt'];adjusted=g['dn']/(g['G']*c/g['C_opt']-g['derivative'])
    for n,y,mask,label in [(t['dn'],t['tau'],t['valid'],'Short transient'),
                          (g['dn'],g['tau'],g['valid'],'Full: export Copt'),
                          (g['dn'],adjusted,g['valid'],'Full: fitted Copt')]:
        axes[0].loglog(n[mask],y[mask]*1e6,label=label)
    axes[0].axvspan(*fit['injection_range'],color='grey',alpha=.15);axes[0].legend(fontsize=8)
    axes[0].set(xlabel='Excess density (cm⁻³)',ylabel='Effective lifetime (µs)',title='Measured curves and fitting window')
    axes[1].semilogx(t['dn'],t['transient_error']);axes[1].axhline(.05,color='k',ls='--')
    axes[1].set(xlabel='Short-flash excess density (cm⁻³)',ylabel='Residual generation / decay rate',ylim=(0,.3),title='Transient validity diagnostic')
    axes[2].plot(fit['scan'],fit['objective']);axes[2].axvline(c,color='k',ls='--')
    axes[2].set(xlabel='Optical constant',ylabel='Mean squared log₁₀ residual',title=f'Fit Copt = {c:.3f}')
    return fig


def plot_mobility():
    """Compare full conductance inversion with a fixed dark-mobility estimate."""
    n=np.geomspace(1e12,1e17,150);N=5e15;W=.02;T=300.
    n0,p0=equilibrium_carriers(N,'n',intrinsic_density(T))
    mn,mp=klaassen_unified_mobility_si(T,N,0,n0+n,p0+n)
    mn0,mp0=klaassen_unified_mobility_si(T,N,0,n0,p0)
    ds=sheet_conductance(n,N,W)-sheet_conductance(0.,N,W)
    simple=ds/(Q*W*(mn0+mp0))
    fig,axes=plt.subplots(1,2,figsize=(10,3.5),constrained_layout=True)
    axes[0].semilogx(n,mn,label='Electrons');axes[0].semilogx(n,mp,label='Holes');axes[0].legend()
    axes[0].set(ylabel='Mobility (cm² V⁻¹ s⁻¹)',title='Klaassen mobility model')
    axes[1].semilogx(n,100*(simple/n-1));axes[1].axhline(0,color='k',ls='--')
    axes[1].set(ylabel='Density error (%)',title='Using fixed dark mobilities')
    for ax in axes: ax.set_xlabel('True excess density (cm⁻³)');ax.grid(alpha=.2)
    return fig


def plot_artefacts():
    """Schematic conductance contamination, not fitted trap/DRM kinetics.

    True QSS lifetime is 1 ms. A saturating extra mobile-charge equivalent
    N_a*dn/(dn+n_c) is added to the inferred density, Eq. (56).
    Different mechanisms can produce similar apparent low-injection tails.
    """
    n=np.geomspace(1e11,1e16,200);tau=.001;g=n/tau
    fig,axes=plt.subplots(1,2,figsize=(10,3.5),constrained_layout=True)
    axes[0].loglog(n,np.full_like(n,tau*1e6),'k--',label='True lifetime')
    for amplitude,label in [(2e13,'Extra conductance A'),(8e13,'Extra conductance B')]:
        apparent=n+amplitude*n/(n+1e12)
        axes[0].loglog(apparent,apparent/g*1e6,label=label)
    axes[0].set(xlabel='Inferred density (cm⁻³)',ylabel='Apparent QSS lifetime (µs)',title='Illustrative low-injection upturn');axes[0].legend(fontsize=8)
    # Diffusion gives an order-of-magnitude lateral reach, not an edge correction.
    tau_us=np.geomspace(10,1e4,100)
    axes[1].loglog(tau_us,np.sqrt(25*tau_us*1e-6)*10)
    axes[1].set(xlabel='Lifetime (µs)',ylabel='Diffusion length (mm)',title='L = √(Dτ), illustrative D = 25 cm²/s')
    for ax in axes: ax.grid(alpha=.2)
    return fig


def plot_smoothing(measurement):
    """Reprocess real data with fixed 5/7/9-point windows to expose sensitivity."""
    fig,ax=plt.subplots(figsize=(7,3.8),constrained_layout=True)
    for window in [5,7,9]:
        r=analyse_measurement(measurement,window=window)
        m=r['valid'];ax.loglog(r['dn'][m],r['tau'][m]*1e6,label=f'{window} points')
    ax.set(xlabel='Excess density (cm⁻³)',ylabel='Lifetime (µs)',title='Derivative-window dependence: short flash')
    ax.legend();ax.grid(alpha=.2)
    return fig


# Klaassen low-field mobility (used by the conductance inversion, Eq. 49).
def klaassen_unified_mobility_si(T, ND, NA, n, p, params=None):
    """
    Klaassen (Philips unified) low-field mobility model for silicon.

    Units:
        T  : K
        ND : donor concentration, cm^-3
        NA : acceptor concentration, cm^-3
        n  : electron concentration, cm^-3
        p  : hole concentration, cm^-3

    Returns:
        mu_e, mu_h in cm^2/(V·s)

    Notes:
        - This follows the silicon parameter set used in DEVSIM's Klaassen.py.
        - It is the low-field bulk mobility. If you need velocity saturation,
          combine it with a separate high-field model.
    """
    # Default silicon parameters adapted from DEVSIM's Klaassen.py
    default_params = dict(
        mu_min_e=52.2,
        mu_max_e=1417.0,
        theta1_e=2.285,
        m_e=1.0,
        f_BH=3.828,
        f_CW=2.459,
        c_D=0.21,
        Nref_D=4.0e20,
        Nref_1_e=9.68e16,
        alpha_1_e=0.68,

        mu_min_h=44.9,
        mu_max_h=470.5,
        theta1_h=2.247,
        m_h=1.258,
        c_A=0.50,
        Nref_A=7.2e20,
        Nref_1_h=2.2e17,
        alpha_1_h=0.719,

        r1=0.7643,
        r2=2.2999,
        r3=6.5502,
        r4=2.3670,
        r5=-0.01552,
        r6=0.6478,

        s1=0.892333,
        s2=0.41372,
        s3=0.19778,
        s4=0.28227,
        s5=0.005978,
        s6=1.80618,
        s7=0.72169,
    )

    if params is None:
        params = default_params
    else:
        merged = default_params.copy()
        merged.update(params)
        params = merged

    # Convert to arrays for scalar/vector support
    T = np.asarray(T, dtype=float)
    ND = np.asarray(ND, dtype=float)
    NA = np.asarray(NA, dtype=float)
    n = np.asarray(n, dtype=float)
    p = np.asarray(p, dtype=float)

    _positive(T, 'temperature')
    for value, name in [(ND,'donors'),(NA,'acceptors'),(n,'electrons'),(p,'holes')]:
        _positive(value,name,zero=True)
    eps = 1e-30
    clip = lambda x: np.maximum(x, eps)

    # Lattice mobility
    mu_L_e = params["mu_max_e"] * (300.0 / T) ** params["theta1_e"]
    mu_L_h = params["mu_max_h"] * (300.0 / T) ** params["theta1_h"]

    # Doping and carrier dependent terms
    mu_e_N = (
        params["mu_max_e"] ** 2
        / (params["mu_max_e"] - params["mu_min_e"])
        * (T / 300.0) ** (3.0 * params["alpha_1_e"] - 1.5)
    )
    mu_h_N = (
        params["mu_max_h"] ** 2
        / (params["mu_max_h"] - params["mu_min_h"])
        * (T / 300.0) ** (3.0 * params["alpha_1_h"] - 1.5)
    )

    mu_e_c = (
        params["mu_min_e"] * params["mu_max_e"]
        / (params["mu_max_e"] - params["mu_min_e"])
        * (300.0 / T) ** 0.5
    )
    mu_h_c = (
        params["mu_min_h"] * params["mu_max_h"]
        / (params["mu_max_h"] - params["mu_min_h"])
        * (300.0 / T) ** 0.5
    )

    PBH_e = 1.36e20 / clip(n + p) * params["m_e"] * (T / 300.0) ** 2
    PBH_h = 1.36e20 / clip(n + p) * params["m_h"] * (T / 300.0) ** 2

    Z_D = 1.0 + 1.0 / (params["c_D"] + (params["Nref_D"] / clip(ND)) ** 2)
    Z_A = 1.0 + 1.0 / (params["c_A"] + (params["Nref_A"] / clip(NA)) ** 2)

    N_D = Z_D * ND
    N_A = Z_A * NA

    N_e_sc = N_D + N_A + p
    N_h_sc = N_A + N_D + n

    PCW_e = 3.97e13 * (1.0 / (Z_D**3 * clip(N_e_sc)) * (T / 300.0) ** 3) ** (2.0 / 3.0)
    PCW_h = 3.97e13 * (1.0 / (Z_A**3 * clip(N_h_sc)) * (T / 300.0) ** 3) ** (2.0 / 3.0)

    Pe = 1.0 / (params["f_CW"] / clip(PCW_e) + params["f_BH"] / clip(PBH_e))
    Ph = 1.0 / (params["f_CW"] / clip(PCW_h) + params["f_BH"] / clip(PBH_h))

    G_Pe = (
        1.0
        - params["s1"]
        / (params["s2"] + (1.0 / params["m_e"] * T / 300.0) ** params["s4"] * Pe) ** params["s3"]
        + params["s5"]
        / ((params["m_e"] * 300.0 / T) ** params["s7"] * Pe) ** params["s6"]
    )
    G_Ph = (
        1.0
        - params["s1"]
        / (params["s2"] + (1.0 / params["m_h"] * T / 300.0) ** params["s4"] * Ph) ** params["s3"]
        + params["s5"]
        / ((params["m_h"] * 300.0 / T) ** params["s7"] * Ph) ** params["s6"]
    )

    F_Pe = (
        params["r1"] * Pe ** params["r6"]
        + params["r2"]
        + params["r3"] * params["m_e"] / params["m_h"]
    ) / (
        Pe ** params["r6"]
        + params["r4"]
        + params["r5"] * params["m_e"] / params["m_h"]
    )

    F_Ph = (
        params["r1"] * Ph ** params["r6"]
        + params["r2"]
        + params["r3"] * params["m_h"] / params["m_e"]
    ) / (
        Ph ** params["r6"]
        + params["r4"]
        + params["r5"] * params["m_h"] / params["m_e"]
    )

    N_e_sc_eff = N_D + G_Pe * N_A + p / clip(F_Pe)
    N_h_sc_eff = N_A + G_Ph * N_D + n / clip(F_Ph)

    mu_e_D_A_h = (
        mu_e_N * N_e_sc / clip(N_e_sc_eff) * (params["Nref_1_e"] / clip(N_e_sc)) ** params["alpha_1_e"]
        + mu_e_c * ((n + p) / clip(N_e_sc_eff))
    )
    mu_h_D_A_e = (
        mu_h_N * N_h_sc / clip(N_h_sc_eff) * (params["Nref_1_h"] / clip(N_h_sc)) ** params["alpha_1_h"]
        + mu_h_c * ((n + p) / clip(N_h_sc_eff))
    )

    mu_e = mu_e_D_A_h * mu_L_e / clip(mu_e_D_A_h + mu_L_e)
    mu_h = mu_h_D_A_e * mu_L_h / clip(mu_h_D_A_e + mu_L_h)

    return mu_e, mu_h


def synthetic_flash(tau_of_dn, G0=1.0e19, flash_s=2e-3, t_end_s=12e-3, n_points=1201,
                    noise=0.01, floor_cm3=2e12, seed=1):
    """Simulate a decaying-flash measurement with a KNOWN lifetime curve.

    Integrates d(dn)/dt = G(t) - dn/tau(dn) for G(t) = G0 exp(-t/flash_s),
    then adds multiplicative noise (relative sigma `noise`) and an additive
    floor (sigma `floor_cm3`, cm^-3) to the density. This mimics the shape of
    a real export without any calibration: dn and G are returned directly.

    tau_of_dn : callable, excess density (cm^-3) -> lifetime (s), true curve.
    Returns dict(time, dn, dn_clean, G, tau_true), all on a uniform grid.
    The true answer is known, so the three analysis modes of Eqs. (52)-(54)
    can be scored against it.
    """
    t = np.linspace(0.0, t_end_s, n_points)

    def rhs(time, y):
        """Carrier balance, Eq. (4), with the injection-dependent lifetime."""
        dn = max(y[0], 1e8)
        return [G0 * np.exp(-time / flash_s) - dn / float(tau_of_dn(dn))]

    sol = solve_ivp(rhs, (0.0, t_end_s), [1e8], t_eval=t, method="LSODA",
                    rtol=1e-9, atol=1e3)
    clean = np.maximum(sol.y[0], 1e8)
    rng = np.random.default_rng(seed)
    noisy = clean * (1 + noise * rng.standard_normal(clean.size)) + floor_cm3 * rng.standard_normal(clean.size)
    G = G0 * np.exp(-t / flash_s)
    tau_true = np.array([float(tau_of_dn(v)) for v in clean])
    return dict(time=t, dn=noisy, dn_clean=clean, G=G, tau_true=tau_true)


def instrument_figure():
    """Return the Sinton instrument schematic as an IPython image.

    The PNG bytes are embedded at the end of this file, so no image file
    is needed.
    """
    from IPython.display import Image
    return Image(data=base64.b64decode(SINTON_SCHEMATIC_PNG_BASE64), format='png', width=633)


# ===========================================================================
# EMBEDDED DATA - nothing below this line is code to read.
# Coil-calibration measurements, two example instrument exports (stored as
# compressed text, checked against SHA256 on load) and the instrument
# schematic (PNG, base64). They keep the notebook free of external files.
# ===========================================================================

CALIBRATION_RESISTANCE = {1: 4.201, 2: 12.4, 3: 48.48, 4: 144.3}

CALIBRATION_AIR = {'2026.10.23': -0.7946, '2026.10.26': -0.798256, '2026.10.28': -0.822348}

CALIBRATION_ROWS = [('2026.10.23', (3,), [-0.165119, -0.165357, -0.164591]), ('2026.10.23', (4,), [-0.574767, -0.574638, -0.575057]), ('2026.10.23', (4, 3), [-0.102152, -0.103149, -0.103745]), ('2026.10.23', (3, 4), [-0.010772, -0.011731, -0.011599]), ('2026.10.23', (2,), [1.59196, 1.58228, 1.58128]), ('2026.10.23', (1,), [5.74745, 5.73308, 5.7595]), ('2026.10.23', (1, 2), [7.42778, 7.54283, 7.51301]), ('2026.10.23', (1, 3), [6.29216, 6.29408, 6.29175]), ('2026.10.23', (1, 4), [6.00158, 6.00036, 5.99444]), ('2026.10.23', (2, 1), [6.57825, 6.57954, 6.57801]), ('2026.10.23', (2, 3), [2.0329, 2.03563, 2.0322]), ('2026.10.23', (2, 4), [1.73297, 1.73135, 1.73086]), ('2026.10.23', (3, 1), [4.91631, 4.91706, 4.91915]), ('2026.10.23', (3, 2), [1.63156, 1.63077, 1.62789]), ('2026.10.23', (4, 1), [4.52782, 4.5253, 4.52259]), ('2026.10.23', (4, 2), [1.22242, 1.22174, 1.21926]), ('2026.10.26', (1, 2, 3), [8.02333, 8.02074, 8.01875]), ('2026.10.26', (1, 2, 4), [7.79717, 7.79318, 7.7906]), ('2026.10.26', (1, 2, 4, 3), [8.05467, 8.04474, 8.03806]), ('2026.10.26', (1, 2, 3, 4), [8.06366, 8.09294, 8.06365]), ('2026.10.26', (1, 3, 2), [7.67201, 7.67097, 7.67194]), ('2026.10.26', (1, 3, 2, 4), [7.75187, 7.74854, 7.74765]), ('2026.10.26', (1, 3, 4), [6.48726, 6.48987, 6.49092]), ('2026.10.26', (1, 3, 4, 2), [7.50668, 7.50612, 7.50678]), ('2026.10.26', (1, 4, 2), [7.39806, 7.39957, 7.39761]), ('2026.10.26', (1, 4, 2, 3), [7.65839, 7.66193, 7.65828]), ('2026.10.26', (4, 3, 2), [1.30656, 1.30855, 1.30571]), ('2026.10.26', (4, 3, 2, 1), [4.41785, 4.4195, 4.41336]), ('2026.10.26', (4, 3, 1), [3.89134, 3.8883, 3.89024]), ('2026.10.26', (4, 3, 1, 2), [4.95246, 4.94268, 4.9423]), ('2026.10.26', (4, 2, 3), [1.62321, 1.62176, 1.6211]), ('2026.10.26', (4, 2, 3, 1), [4.72573, 4.73153, 4.72933]), ('2026.10.26', (4, 2, 1), [5.31577, 5.31473, 5.31083]), ('2026.10.26', (4, 2, 1, 3), [5.58178, 5.58112, 5.58113]), ('2026.10.26', (4, 1, 3), [4.99068, 4.9923, 4.99269]), ('2026.10.26', (4, 1, 3, 2), [6.03235, 6.03312, 6.03394]), ('2026.10.26', (4, 1, 2), [5.97412, 5.9732, 5.97213]), ('2026.10.26', (4, 1, 2, 3), [6.2392, 6.23851, 6.23819]), ('2026.10.28', (1, 4, 3), [6.58537, 6.58395, 6.57533]), ('2026.10.28', (1, 4, 3, 2), [7.61353, 7.61507, 7.6149]), ('2026.10.28', (2, 1, 3), [7.14965, 7.15151, 7.14592]), ('2026.10.28', (2, 1, 3, 4), [7.236, 7.23402, 7.2322]), ('2026.10.28', (2, 1, 4), [6.9016, 6.90128, 6.90119]), ('2026.10.28', (2, 1, 4, 3), [7.16612, 7.1632, 7.1648]), ('2026.10.28', (2, 3, 1), [6.07756, 6.08074, 6.07451]), ('2026.10.28', (2, 3, 1, 4), [6.16605, 6.16257, 6.16485]), ('2026.10.28', (2, 3, 4), [2.22144, 2.22096, 2.22077]), ('2026.10.28', (2, 3, 4, 1), [5.38386, 5.38814, 5.38781]), ('2026.10.28', (2, 4, 1), [5.79555, 5.80011, 5.79469]), ('2026.10.28', (2, 4, 1, 3), [6.0748, 6.07371, 6.0738]), ('2026.10.28', (2, 4, 3), [2.14298, 2.14188, 2.14034]), ('2026.10.28', (2, 4, 3, 1), [5.29583, 5.29778, 5.29423]), ('2026.10.28', (3, 1, 2), [6.38247, 6.388, 6.3868]), ('2026.10.28', (3, 1, 2, 4), [6.47475, 6.47481, 6.47215]), ('2026.10.28', (3, 1, 4), [5.14978, 5.14953, 5.15039]), ('2026.10.28', (3, 1, 4, 2), [6.20838, 6.20575, 6.20476]), ('2026.10.28', (3, 2, 1), [5.61182, 5.60693, 5.61604]), ('2026.10.28', (3, 2, 1, 4), [5.70024, 5.70303, 5.70018]), ('2026.10.28', (3, 2, 4), [1.78556, 1.7886, 1.78139]), ('2026.10.28', (3, 2, 4, 1), [4.91186, 4.91204, 4.91344]), ('2026.10.28', (3, 4, 1), [4.05775, 4.0633, 4.06286]), ('2026.10.28', (3, 4, 1, 2), [5.13796, 5.13854, 5.14024]), ('2026.10.28', (3, 4, 2), [1.42927, 1.42697, 1.42911]), ('2026.10.28', (3, 4, 2, 1), [4.5623, 4.57078, 4.57632])]

SOURCE_SHA256 = {'3.txt': '96d925f6ee656489204d37e4fa02933aa43030f1a19723e891d301100f2963d7', '4.txt': '5e4d2fcfc4eb6ad3f0cb3b03c0f273832727bff17bc2373d2fae98b47c718ac9'}

EMBEDDED_EXPORTS = {'3.txt': 'c$|E_Td!ovai#ej1o$8HNMVq{?RWQnqbN3@fuw0P3C&w43Z)cT5Ls1SmDG&>eXno1``&R*#F3~XvrZn5aPy0~U6w8TyKlaJ`=`(U@b|Z0{^IxF{N=|#f0=*&?GJzX%a`xJ{;WU#;ZNUv`|sa>`T6J1zxwtse|tObNBre4{_dNfzkL3mUw;1X=da)W$9G@<^YgEM_|soL@%Q1Ezxewfe*EtHKl1CJ(%;_leYBd=$UXL2Yy92+70UT#fB4TIzudn6^i}`!%U}FIKYso0+i(8-`L{oO|MS;xzW;g)_o!d}$M62_%h&Jz^5ygIfBf?8=YRU`-|tT@^P~RnA3pk<-+cegpZ`;H{QOTp{Ogz7AAb7g`=7s)RKNVi|NQcApMU%1&wu{>n;*aY^!2Yl>9fCl{i;dd!r%PrAHVtOk6(WN`nQtuyKjE}n}7W7kAM1lJ2d4V{_X4k_l<t{U%veG^LIadfBS!b`-e||`2FAfN4*~pv;5Uxf7U9DBi~z$A*P=7zgos0zDl&;|LMz@uWRAH`TpB4UbkQW{J*tK<sOfbM#{0(qvv1#;-9{J^RuS*udmnnxa%uIZ@u<XVyMjv@y*xIvHbS)uYU9EZ~qcCje85R=aJHpM+xb7pMpO4cRzpq{C9u-^S^(7JLgZe^y8k7bc7a4O8Ue<eEEwO=o_!luYUWtx7h9}r<jiZAh5pk$A9_DkMl!PEc%RA!f~Y0jw&&}`~GojK8;Z&U<$E}GVVi5^H1OW|Bt^BwGgc$QtmaEti7yrFFyZ|AHMtk>*woNweH&0A3y*77ws&obN=|ZKm7F5mv6uR7pqs&s9i@h>R{%~>V5h3|Nj2F{}t!z{Vy$F*KhSYLQAprQJR)dGYOx5{_FRj|3T;XJAM4W{P6APU;PW8*s?Yvm6CPZwUxj8#r~E5{@*|UO3V1Ij-Eb~`%tuXLA#uCKlhpY_W2K5P@Rmw{CErT*GK+!JYxNukpElB|Lx#ree=`5>l^?4^_xF_@wa}`bu4}4)^lZR|KQ!!KKzeA{qXgNZ-4myU;p~;JZQiF?vJ_+9;DO%<)<(F@bRO5JAZjC)bGDzLi(WJ{P?5JJgffAPd|P4<)`;+{fBS<>M#6<@Bc&BLN{9f{+qx4QNMU?!LNS)6#n+6gYd&ojQrq}-mmEX_uv2Zd+8pn`4;PsU%uy1>I$9*U7L6Q`TOtIC;svmS{<FkI(pZs8ofmw>(J0)M=1K25^G8P>lnJJ+LTgj2)VyJgSUVF=AZq~%fHnAd+U2kwbazc`Rn`>LegLAZk2rJe5|SEag?gh%cIwjYW7@?Bc>`PFhb>Hj~+GYJkSei^ueXpiEU9c)vtv0?MS2v97Vq$M?S*sQ?B(M3WqF9rsJcLI=3U2qqR_T_77SdqqVNZ%B2N9JRBNXj}k_$$C54R3Rf+o&gfmU4@)L}c+8S16uvB{u4`0s(Dg|p{O(iX<J!?mPW2cX8wwLEY{)SXU6NlirEtWoRU0Fxkfc3|)oRxl1uba%XlRq(<)e?5M~v@<>h)PDEm+Zgmqs|G7Y;pD8FeW!gEVYm<UOV^mY7X{jQVzI;oOdr+AZCsKvRz)^(AIyrY-ASX*OQ}#K0pS`otRQC8mytRH8JJ?z>hm^l;ZC^jb<2Yw<oBDa91Sq2JI+9$M*lV&?WF=8&4d>ogsq>(iu*>lQKv%`vrnPkN`9kX;vqQ*+e9XNI}jRo{`)l2J>eZ`LXg9nz#v(x>+_deSU)?;eh)&UI`u15(q)$+vvyT-2($Y4@_f=jBKp*FJjDA#Onz_-QNC^H%Dd<a;=(&P3MtNj+{U5ACehImdg`te2ES%h(Z3I#lPVm7AV9{YhGXjNg6o^ou@9msNKV9gmnt)J01^b&t33=xt~-bSqHs+RGaIJ!zpi$W5R9!APBI-CNB)<s1i>uJvA|Wx9?&ujv%?(DTb7ZCxLC`^4U8>9tKdVmf|XOr1yVaMLNO`@|y9V3l<mwZK~NeD6om4L|x4lGWr0_AG0mhuhGp(|%|Dw{};@T2|E(viExALqF%}rdw@|n=o{YLmP2?G_>jF>iDFpWsE+rhqR%j>yot9ZXvZ{T6;}W8&dU2S;u75qzzzh60sm+lzovj4PBJHMU3c#d0BMKq@<;XIx8`n&eIZ7Mwu2~Um^P>HfblmO^Ym9WZcwK$Z%wxYm{Tvj+a}J0_b{crBzs;4@WknmQ}{K$s*8CB!|AT4IOK(u}x@QYi@&7e%5u9?!a&evJ9j@bO7c)>RIa23pLlfc4jHqkhI`v9^IP4$t@C6UoSnROR|JaGFLc9P1`TE*|gJ54?*ZUVk2A|sfQpny=0^5yIZ^Irrfm`k{Hb?#reU=rkkW2r<<(3JNU{hi-s}juw7*#wSk(LmRTxGk808rsT-w%+9v(Z7P9GK%zB7*CrepEYR^#!lChy3*+(6<!t1p~lxVy*GS9NUuiYhvUMuYFLki*&Q6{Mlbl2ZX5lHVuJ;8@$%qR_=@()Jp&rP!EW!5$oG%5F0i!2%WMqfxBG(8tm6f&8)Qb{UJFRgwq)tQfTAN5mytOrJie~DO{76Emj#i;sMAJ(kNN{zKGq?DyB)rcY(<yJK-J?9#DbRG6^2`K@a9@wrA<sr}>-gWVs1lG0EP5WS^w1<v}HbfUPM?{ig9&~(o^f)}8gbr`vQRIP5!F%BsQ%|I<&$NY1yz$y0+5g9SK{#0nNp$kHfqFOW)N!ly=PjhPpe7=#Ne3fD%~ehWFTo*oKJBX|rG7=1NP4bItxL+%c(LWINUHYnNZoT@Y;1^ZdnRg8Wg|&*6f9>?B3GGMox`LDi35j4ts|d#VfO<fi$z*bug}3F8kdk=S}h}_HYw^p)gs**b)UL)$PzMY;yNihb9w}%x{LN&M~kU(KQ0-i9<>@aFSW)u-F8`5SqDMqU+ezCz^eCmjDcd&k>MwLrgS3pby9$F-59NT(QVR=l4T(gIha!adO<?p=8lf37D~%i?z*k}KJvz+%H(e9L6SD+tV?0bZj5E0Rur19EDyXCkY=iDqJ{2yOQfyiC1aCOq)(Jk(hrg`w7c$>98}qHTK$hkYU_1{_5Y}U%cyO~qAy&UY3(GbP~9K2^uUtRS=GiT-c{Z5C1OY4X)*Llg>8&Sy;N*s*N4W&!^m6dcU-!8OQ-`Vizl|2b?<{Fs%^r_mu;p+OP7dRbu2Zu09~Qy{_a0ZUJ1}X8YuyJ2?v`Uz2e##=@%Wz%=XR;sm;<Gq{$>!93cz(EXlCN(%w>Hx9X+W<09p}q{MWdhN2FcRGXd>J!>^b-m3jZ>Zz4IQh2fSfJZC^%|HgCZo8JLU1AOt6N{xgAgh1y*y*Sr9O6_y8d$kTsD#9ejgZnXT(u@$SnP2m9P}(DqyLq<ei_O@MpUu>BcoQEvt-0pp5}4Yvx%XRtk-l->baA)K@J=ttwh&M*&OM&%KRV906jrnlh@mg?J{n$vnQBMBgzSgAC1KN)bpa_s5P=drvuKa>J;maw`^lXJF4@mBXQU?2)rKwT_xjIr*unMi=Ou&=TSGdEh#gWg0-6@#~zfFIy6%FvdgxVQpZyIIxw2QRA1&1(_5h<bRg0trS3t)Z`2u;zQ5@P+@&KDS5!uH{%EAsfo{EwL0MwOzO>#^9V}#$5~yKHDO*TK6gx<_kbBjWZ3?-B`p=fK$fnfZ$-ho%nRI$9WUui>wdw?E2|2D>dS0nh?}UUqphKsVaB%9hGLi-Zx@dVXDYaysm#yf?CMo$U(+OqEM#}d6Xk_g~Fr@?oXNKz|8&c;g>3QSK`qQ;P{1bh{?IL8HH!TdFZwXsU9ToX&Nee19wWO@o*03}U?yih)DQEqQ&ieICMl3;lPnuxB-_&-=!__s3D5M7|i6Y92Dr|{vW$`uCG34p_Xk48}5|cVA=St=lcC)k(XLcUfpr=9VRhBsBES5GpPD-ci?K>?2t)#wBQekax^@<f~8TDupY6v_)5=VzXpT1-anAh6bfikvTEYDIGDB+RoTQ3<CtExFiEGFH`EN`+z9DeB*xmh2KtR3TDC-cB%vVO%3kcGvXM%#G@ADY=i=|$PIm@Kk=q@VOM1U*Wpn6bz&9tl%n>mH)Ah=O4uvqWn7#@|=HZ2HhGrhHdDAb5Pb=>xN|=^Q4Mp3cH$ymFZ2a_D&Luvtt!BRaf952d-{M<X#p<%H@PktI5<IGt&2X68~#{lqC}=dvTLu0QXHJk!EOl(nyic+06}$wjwScRMWYqxVc6IdT?h<%mh|gOv`uTU$=4R_!}>sf-S0oADbu7lfZ*a_Z{n1tMZn^d#Q&F5a<GbT9Qb@_xP^sRs;4k(iDC-e+`Rjzt{Y5ayx9|E!Lb9K*q;h~F5oz@x0@EoH#zl=-abl(n1s_T|Cpbm&c~2`O>x3Gr=DL`JUm++Lv6Aitptqt}>Vfn<L;v;{5fD8j)NPm)RxoKC4!WNjY})LpB1IvFpydgpo}j3?hPNjJ~^m2I1eD9V6CaQwsM=UHE<lW<BYgOV_ojJIQH9}bR%)J=%ks63Yq^`_50Wz=Dnve$#4O=3S}N9vSB!j!rzmyBBJ>T8F_mx_%;Qt6SvSC&os(Ma}CS6>^Y%VLeMo6IAk#~u$buf3G#sg-X25>6~p133*c35n%zDGyGLKIbkgbMIxscEwZeoEE)XT}GHzImrLH_EIL3%)LNvMfclb<G@;LW#tw>?4>pxeRxO-={=F0<*wUVZCc6?{%E8h(T4LZ!@?mu{j$~}@(%F~%@<SmS!W94HxZ<eoCR}VdObZnm|Q{~Mmcv{fF-7Wuh=2ire@!FQgG7rIXl`n#nd?|x{(#FF%Yx!R2_t5^{0o6*$PwWL@A^IC8kUiyiFNlk}Q8Tk{BVDrn4c<2Xwk9GRIKu^G7;gJlvTZEVtSof8oHgdGZaGoH_;B8cpxTl2Y$w%nmz9F6}@jB9aMg$)8K7l(Ma*-in+Oy)>GuEW$*T5oh4hF6f-eLByllGYV2ao|RODw8lpxaauSs(j9{kFjDH#)o1B6%Sr6>Inh4rmj(v8G!Ihhd`g)}`)hvVlrj=c35Ac{a$??c8Ds~@^yS!+e!%6$72R?sV)!ye^v!y;u^psHFp@g6f7nggB83cr3=q5$nL~GJtWe2rMX%7yk%W1)PmTagrK&Zk*gHCdgz~1C1({}m8Z~?HCJ|(jDma}I>2STZ(hoXDaVsF5EgxBVFj~gw1WeV#taq=T8YtmubQIVi)DoqXr9rc^BgdmDYg2>_&EW_}3<0MO9U=R%=8s0!+6jl~+DnAs@J=C@qNO3_pp~0a>Rb{U;=$lz%0-ch5v?^&BTtV%f#6|BV?UYrZ(R5!y(VF#g-DZR2Tsb4+%Zk!FK7ZhDqK+VC@4_v$b*#4aSy#adRp~K#~=zQC8BlH_kT2!dzExunu8JEII=7H(n3I8%QB_xxDAcxNh1(bODUn{ouCO$G5s#}E^jBSXoy_R!@DUZFD=ni;bGHKcg*T`iW%kX>&nR(tb~+w8Kws-6`HQ+gp`LYNuq?ozj#!6^itY7r802)M<WTTX>D}B^eXxm%5c{q(5CD3)n(@u+iQCEyRDf(EK2Xm;FFGcZe}A7qbCS6U|r0N*}+oD#caA01n3jtyT&`b-LHTfq<QhECr2^m|7aOxS>vU~2QlTr$zwbSjL3u0+i=$_r#W-Et)O3y)F(I08iH2(OA9qfR(G;<GBoCobTXwa2=428IRz1Lzod90A=Xn)uVASdi|f@Svt%<_$D4G~uxXvN=ES&9J5L|iaJ%*Jjw8(MA8Cfn1dr?0lzNhe;Hr!)v*t$bFmw7N4`r<%jl}ZOjgW?uo7wzoYJ92~zq01@SosSXwA`Rf9F-)*0IPtvGW1t<rVgYPSPY$c{LFPZxz<Q2TUg&kBuKV{Hg7mZv!v8bCACA`U7tt%Oj{q&X0niLe~_}EXi3;=J@t!+KAFJcL0YA!k49FkHN78l5w%%$y&X;h<sy$NPd)1wS?*VKxJ*(yqQvgVnGUCrq%yUHdM&m&t4AU^$08$+ua=3AnbJV0$Wb}PlyBS!Wk_!{Hj_L75nyQrS@@TjvM@PM)(lY^fAl+79_?D^JM?lSX~<lxm&(+sa^Ojpt)-npW}F6sY<P|m5~op@S;ldy*G@)UG}(FDwG}bQ<X^0RNPy|vbC}N>D4B2?KnJLf!-2KLlPnWi9*C|Dhod~ms11-i-KBA4lu5>_&x(9UanksMk*vO!TR#d!tY#D>01qMY*y}p4B#9F7WZZEHsX59wmjzMCA#EuO@mcNqU8ib2nK~cEAFGp4<h4s1I*!*1-um9k0~0h7{MJ<T-9!qd7YUCM-0KfVme@eUSC;l`$XPTQ^7MKxdiZFd&MFa(gCt;ZNLFu`HmqYHYITpbnUz3T)+^sIaEF{Zxs)W3*v@jw*>K}>yOS$ka_WdU465@e!y-DC%EKH8F`i;({zk7+BPWa<DV-TWC+I4b2Q{0lsWI`Tt|Vjz$^lv;hXt*8F%B)GZIb0AwdfcIkr<4E)X1KNbu9^dbLc=66*f~3w2XDhsg-P}ki5M*@x1&ix-OI5Z;1}GmQiLAGFID>zR*iZEK;5S?knZV>g6Aj9gUX9gOsvZq;I6kWC7^$mf=k9;E>@$Lgu59GF^26bu*+h)7(p~0?HKSQ|9*S&!7i}wOQa|rcTI1M-4yWI<o=INqSb|FEP;*EMIh3KvpM#egepH1NX%|u)j!ARh*=Qlcz(|`3H~UAyXo0u=TXoojHNIlb8ny4Jr2z2UfhhVB;Z7g_;QKCCtl?xih64I8X9oG*EMv538#|;vLsDZz1K=>bXAHxxO9}qJ|DW6O%;<DX%PbPf#$No(YZ-_mbP&<@GXIIg4YEZ<N9Epk#Bxkiz_zCX-Y$wltN7@;hb%@#RQKG)Ptx{Yuj^(9A$ukQ9~&(U+WtITQk%on$*1Br+D-yWl+Z6jZY7rjWzfE|@e_w^PoWJog_BBm<!k4LSvlWEdmfs}`D&A?}h?Ml|f@K?7y-JtTM=IaUj*gA9H~^U3+6k#aw}-~EW|ICDVqn`A+g7XW54tp#Zu=_V~A!TejlV=>`mli#kV3X(G^Egd_GV5wRVF23RtBJkolBB(3jH1!S~;H1i<Ax#dTVM1k!oE?@%d-!lYWp0Cenov(|AeST%T+t?!i9~rdvLfE$z>|DYOBL}bB47>7z9}c!U5>VXP>~bgE-_R2DE-tzQVu%2z35zz^>E4(AHJR08DRP71W|)_Lt3bvV(N;Jw(44AG=2b{(ZB{^YKI3gWg29RaFmQ(%K-8*qWrYn6S?(rWF?ND0h81|{S#+WWymESv6-!vn4M^dHkJH8Vq)BDY0-1x+D_s)(hheS4*P<3UQgr%vw;#Mr-~bf{w}ASvLmY2R-0cmfZ39=WLp}U)dxANQ?Z%Xi385JRl+SOuv$MFNl+S>ny6Rx!jr}UG$C7s{Ngm-D~<*Z3XmUT>GUzW$4+vUu?0<Jr1ai^#9O!~eN}KMn!p{Pg@Zg8$#Hy<OHiVLxW|KFzoBMML+BeMeQOVbc1#O$+6Mtc1a<03msGkve>74<DlrF4#TmuZ$DmlVRIm)GvgFi<kQdgml9~fEB_A=82Gz=4a#oZ;(YsN$Ra@C(D<T;!ASH(JGzf`vyydLJ8F<-l@(pd#BqEerYyi#W;ow#_J{nSs_9pOM!A53)PZMACYNXC|ow2;6=hBKu%P<_PMI|ivOjnnypUDwN1l5Cnbo|7I|8feN9m<NbyT&@XnnZLoT25QW)W?AWkeAt(ehR8Loec-jAbmk37!y$+Nh-BH3R<*2ylFr{P(+8E0b+xsWcr6A^)3>~CIv&L$ATt4NKcv!rIwb))uBD?5O>sEze|nV9Wwypd|DDp%;_8X2usk|?Q=$~M16Mvs|F}cTNuxSi#>>@Pb5p7Z_~E}lpwqG=HOPIe-eI<eOOLOML>^iTKs5a<&5cU4;1B$?&&J&x{%k9_+iN@Uz5D)Otv**O*k@@)b^!`oPU;b0$I51A@IPCNCzH-Ko)>B>%ep8SZz7+t@RocCM(vPGMEaHAKatyD5pHVO2|adL9SEI>?ahsNdIW0j<zmU20@`uu$USOJVO624YtMfFi^X)qHmyl{DuM3i>KL#vniD<<5wDMiK*RoL>kBin~!5i<cM_mX(wRBXrt^!kzgR=ts6B6l<SRquoD`o>H~B~Bc&ah{37Y-3RcUm+{=OFTr26k2q<7liS~<t8DpROtCNiZ<Z!1ZE-8rz$(ifi(DU^K7^>4Wh&^pBq|@SvbRdpX>deQ329ox_JQckr8N`b7Kt{{TC?@Bg07O0*w~Zj3{06ynI%f!(uxm0PQ3G*b4&*HI_DeaXHeV2%8ih0knf&#7mWkB@SSuY=?*tPuy|tlIj;EMfN36WNUi&R3afJy>!Kjl}VGxY8E)j!E%*JLC=$9@vT2X2+VS(98dx-B*fRV|?4Jl{IDaAjr1iBB#0KObpv8A(@LjLrGIS)W6r1Wq;^O6()RqCc|^CV^`*#N>$N3@(*GlB-t>v*R~X^FX9bl@9-2n0zKc`3w$PLQQ6!vPeGluWY8N|Z4<Qcn=~K}Z5A1(X3*2a*!ZAt1YgU)(+zDE$Bsx)S{y3x_{z7@9*uThCS^HoBhL0Px17(sVaCha`VbDal}KM?u4_&E%LnA`SdbM@Q}o4@w4#8{+Aaf#@N12f{+@?eL-IHT1<&RF94)3V)ytb$oe0TICJY-NB<s<D-$fx|jx7Hg%fopp+`0q-^Nfl9No2P8>cLh(?gm!0mOmbx=<^ql0wgE@NQ1Dl$~tm@p$^Q>3U#Lz2EI#8S7Mkq5sZC{Vc@AapIW0iAUoGVJ0SW#y3Ep9fZ_KVW<Y0+=txRNKmlSmLPD1|amL-y;u;akr$CqS2kw{lcu)!gW>&WT&q9Y5(XfkmWjv$FC3x@rDW5s3aR?r;G_d3?L&;1uL(YJDSOj>%iZn00y8hVq|?Tu@f->ZqX`%cnIYHh1@w=OI@MoK@w7mW!Sk7GZP&E8xh@cBy$IANiL(1$t)vY0H_L#AYNxQfYsAM1t`28N$HU>aT+|Y0NDUsMa-&RPJKfFK@MtiB$3Fk6goczYP0SMZ><y|H5yuc*-;kG$qz8>{l%EdccY4A?~Q4v2HDpKd&?EsFC!tcc=+{AfER&hqGJiGUcmw!kpmN5?#?oFS{{}RcaNyJu}xB&@?=jl$AHH!vnrCGOgDg$xf(@yten*0Bc`K`XI8y$zza1&U8g~K;du}tLneU3Iupi**8&#|6b~$`yiHf6P$5OQ2K|-D7J*GUC2OR#x=bJ0)=M(dP*yD(O%yV4<z!De*L?{GCX6J1?wm=$kO2sS>{2xLl#Lh~xN8Fa7VU8ITBk>={9;Ts$dxy?FYg~{sTMK;=IvphNIk4_MF6){@VEwLVH1FUzhsk*DSu1mPT$Boumz4V6E%#+B7l7X;;74d+HF8)crtVt^v>R>1I6qQ9VqX*Zv#+AVf~`b0Ub#l&W;_b*YwpOo;`ypbT0G<OGgx^kaP(0h?hE<d09EqL=_|*H;;NKP#5VmzJaT0>%yHd>yoZhkArm5SrnIqiDpg&9ZzkNS^Wt9glwNqCC)9?3zeu`Skht15!WE6W8r|B6B+3B#Rz#u4hJ^}rQ#G#ibm!or7LwuIEcVeErSEoG&je>LUe@oIA8U^>!BC)PM&Lt23pU$M8=z*Uwx)G?&u`6&Ue7}yu+w?i~Yv1a2yk}Qu2}w=q<V60S3gv@rGkx>c8WQK~%6L?5XI=7K?{=iCjAjrsW^%JV=s+t5GmJS_X}RX+V1I6i;3x)gd5;*F&H;$JyI}B}NVtRf)0BVqMblED$Thq0$RQ=BVSzb$s-_bJB0<t_-G3FMi{>B>%#zEZz7#NX``T7&+`lXH4mcU`5hZz<ZZ;QaQ$X5{AKa)GCwo#KddY=GiQlbh@oeqfi*LpqGpQosX&w#D8Nx<5TP$?x^z^T6fp~u`vMV#}=>6$cteVh!=tD>SB=p(<}Zwr~-S6xSZUSdF^5YYFS)lSeC(wZG$u;Cu8-sk5oWOr|Vc;*G_s`mc6XNvS2&lVI4L`aw;bimRtkxMUH$+Cq<vBT52|IB^|+1poFdj+)9O9Th+yNXw=e?qazwjnB~PFz69@%6h|IQIueR{oEtH+_0SC$tqiXd78vg-!8+squUDK9C%7*l@Jl>-n0bOq@&9Gz;&jNiLnmL-5#J=1jiE3ko$JfMZs;OUP0s#mm1Q+KtP(68#Z**;@&djXq~GNY0fq>;*U|-Y1jJjxUoNMJ>~6I0IO@>+3VDCImfgPbmhPbbKp%dWrrJl7QW@(NDl5rbS}V!ilL|b&1?qxj?sMZorQHDB0%DMM;E`X_1r)wZ8SN-YGfVda>(&;&7*gY`sM0Yc*)U}z^eZb7gKt5TJ%UO#o#B9N0QBUsf{@Cmaf#^x-pc^VB^t3LOU5EM0-Tip%pB8s=_t_chF>S1=^Lgy54z!oX>Xj4AQ{Cp>sZ)0w2Kd!eD@cFfL$cNP+E|gO4K|^^+xG9(tb-g5M$<AX!d41>d!#?N8^&p8!A~bS@*`d-UsmxV}M(lIcr>jD?*y^If`v5tf(7tA#lOX7}H}(C*v#l>Jylq$FgFjXfpwAJ`ci=kewn|4jXBe>avPpA)>6Jr7;P?P^MF;pQT!j96-mIlW_f#?f@Gi8>{cSB<)6b0GN7EtYc@3hP+R}k*DoU1R)SZ(<PO*yH)a&`g5Q>W$6@-TPHDqdekxFU``}G)kGZTr$KtXq|qp>EA0FfPrk{aO*#mK&s#@Ys1lKB)Ho1UqILj#5`5xyG%L!NlHss3Z3C;>It6{NyIO{!U7&HsS$B0%JV${LTMZbaQsO&^@#;r;U%N_P1B<22_s#{m!&No1H(m@fxKeL94`rN27iq;~hEHQW#fwzaa5%ElG1SVC?<e39&IjQD_m+9@(njTuV35fHIf3(xQem#bK!>#JJ$IcfY3qs$(r$$m2!KRf=otqArU$=3MNdd7%cm!pJ5vY|KaX+BfwUgp&jeQMrF8BF0FN`XQ$D#)K+`FAar%&(hZ|9@2Z58J`GKiCmP<aU^RRH-D~h9S>J+PzdMOrCI8s_ct407FTOIUhG~Rn{UZA3kw1#tlH>pQ#E$6$jJYNjbQ;wF_eZsw%9|I&s?{oti*;X`=O9u^u4hI1Qy=b!N38w)Y%;}O(o2ym0gLYf;0SWP;AC5Cql8CK|Wt%T+MenZWWwVzE@{NFVg@&SMmzgj5T%5owAotbs@r2qDlbU)lsA4|pFcY!~6NCzUpgCN{VUJ5Zjt)nbfF*D#a_@w1<Fxgcct9l`Ol})xOyYDpjk?d7++rtotv8;i|5-Ip1r*jpEKb<T523IyKn2tYgf96=Xt_||)qu#SPo(}U5Iue|NIDbqUrVe1nbsJ{dtUIFFdJJyFyx2pU<UY_MqtYJqxSULzsSLpV+jO(38-V~aIWzPx*f?idsNivCAgF)XbkyCz&-iu8~!mqOl8W;n%3z6OON6sb#>z=Xi3p{`eb?I&w~JB5U-<-W){k0Wl#ug-Im9`eUxt>uDa054FYGnYoyB`KvG1ssXqpSy_HMDF|l|uUo*8uuJWq|Q2vD<0-%FG=7*j87aSmw)sLaJ_m9H>a|$dR&EgS|!PSEg<i#LrhFCd35Mu=owY`H)4m{;6DJVk=y$|5Oxfs25t~I)mQm#`#(4aDLh^(!bedyOFdrAVzAdlLW2RaXWdB$mz>B9qcBCZ;dJq>tGI+7w*az=amKsPwE7)cyTzCL^{WxZ%Wp9kq1G0|m%*@m3*Vb+6NZ{(qtc$6H;nkQmG%#RevAV!seZ=T}mE)bKG?q8UM;PO8oi;k5y0YM3^4gbo8zKCF#VI`q(LbqNDk?In?l_AWaq4!Tl5VOaysFn{ajx1l$L3<vfH;;q?&#L(cX4n8%FOW9TnOu(_R+^p=vhduV>H=gU{W58<-J{$~e{nBkc|(}j3(~0#XGQRY6IHIw+)9^tRN_Qp5SoXcM@qR2maO%di%fEul$F+6i%N+lZWl=a5DxP7Uku_+r+%uji}RTRw?$P+bH#t$`^VWJKqQAcEEGYYrmLZTmw*)mU8bo_n-rIPWcj`NRP$l$xD1z)oD7}UDIcj3?Fv8vl=cJ+Qfoobe_uaLH$8ZmpA=V0zJoZ+Fa~<O7(_O^026T;;9MqB)N%^#;Lekk(M9+|sJ85iWeO5Ojcct#yTk))He^}KAX@TyLY9wkG7i*CZJIBlUGfnC%J|&u^9`_TL^J8Uw{3YwryfBgX`f_y)@r;foIjcUd>#ZiyyKb!d0rNt(YA$35Lq5!YZs|m=UNHqyRb{Im(PlRHJ|cPIHy&<hkg0VM3gt2gh<&2vf<BbR@c=_4khBv0ZGIk>)cH=&9P;p00wg^QHUnlDAXhr2x9!@Akg-NR03*#nul<W==n9#MLo`|N4*br7SA`3qC|$J@XY4aFWIzH2Q@vlE#YchJS@#oM4;f9)ml){ao5&y>#`1RbrJWCYF+>kS`uRirhJunKa?*z?p!^WnuGrUdtP1)f^&g#M|nvtPa6{i(P-46OM`~<>S^A5w38aQPv9kOxfPqFr+ly<$YdiOz2rNg)fHn0NB~JM-l3#_uI&p(k(2#3frgnb$dm)9O{aW`LbAi09o74%J?*6Dc=28gqUhA5B~-bt8%&<dX%3m-E$#h7=R@s-P>?na%HL;#t<~#rYWgI`?qr{qge0J6kq)^hAci2hjQ>LKT@nV6q)_%q0jHbHDPA?e+my$gClnAod?t8o@-^TKBQM8;9bOCqKi8<*BiOl~7<@CCqvQel3Uy>KOom5NA69OkWSrs0kUXcB@lr3Qaq8hJSPWndrwNN!ohYa$)af)17w(jhQcLEdC5OHO)08?!9FDZT2f*Bh;F6lRdZEUdf^}WI{$dQ!b};kPR$5n(q8CCLq#eSRk8rc2rURxXlK>GG1jw_iPKi-)Ng%&zGv>m!to;K`kE7uTd_|?7_s`Ma;^2f_G8s3GVqUCmV;q-k0jnEy7>hS8GfG2(Z^d{qs8apQso`tqrxd~~2bIDOUo7E>b?BsZtOWB=p=~8`an|LM4wz1Z|0Rwk9g5MPtdjzYTK9{$570_zj{;`w3wbb90_IzyR)TzCm7RFx{W!q8AO#&K(izhQ);03zE@XWgRKYhBA|P<FEZD*Zkj5C`=-1=NaRv%cf;rzl#ruYMQb#{72@4_ZnC@{=Y;X$0)WSMZmeEiu%H9g!ehNr@21BY7;L~Ol7bH;j4`OETpe$BcYslv$3rHbu#4czcG@k}ldz>U(Lt8M>qnJ1V$n!`Wt`6=G?1GLZj35XIR|RIEi~27CBL!Eod5^esC}piPvQS8}QR$)+-tGiPNa5-N_9Y^DB0yA;=i|oHLyj$D^De(MQ?9R%foqU3*nHrFp+^B#{XB@0SBm$Di^av~SO|<lCsQa-NJ~Vh49L_HZ?T8~jHxv(yj-V<BymOvcM#r7L}<%qbsX6l&^^X?;OJ49er{rOwc((o5?rxqW{f-MeL2o&q5zH{0|denUPkOFS3SKNLjfP5<HL^2+QZmkGw}?{I0fv6<|j-!uHm5UpA0v#i%ZF5qI<)%cI%xF5MEh6<5ayp;hs2H*%iiHK3K^JC6Fr7EoA+;Lcm$WKBQndqE8snVq2hC*D!Ax7XCD-n!5<|2LEJ!%q*8cC_KQeL<>L!t`5+*Ph6t`TNU{?F#V^1gy{-wZ*gBk%Gala$+0ht<3?m5W6|mYFCx%Nx+L4B<an8D69^bc(RkX}dZy&PumK?J+Tlz%*1Nerp9fVaoCs}T3iJiRz?8|HJ~&$T7E+WIq0^w8S`rpY(jj}K$TXc2PDs1%Mp!W>T{VLKmQdydF}}(fFM#c_mo5oygH=;D@T75;G8Hf=vST~BJa5u59K7U=AINV^;abqpPlM#xC0jD&qL-IMs6I`gqS+v)gw^%WY?=%yUP8Hh6dSp)vs@CA=MH+u$CfanlBS{|PCZyp{K#f5$Z+i;IXIbyvVS2NiDYIQ>Xw@5^%gp@29?=>>|^2xI>`t^J0NMF232FaU6zzrD4?^_l#!7b+R|wo61fc>AEYw}7$D|Oq3tCdl_aQdec6%@h}%@gaB_e=b)gm1M353+(xI%i@0wYIQlUeoqD5t=#~DQHfDv)I10Yd}2VoR`Cx+#VL6vYtqK?0oGZ0Sp0ho&8A^Q$CJWxo%dbg4Rs!48zNOlgVgkS`8)y$x^B&6hIeh93E6ka;zOa8D@E3^e<Ws!||YlC75{u<)7i;KnrMx&tu_?S4cWsnWR$A$^<JP0Ukc7ltX^m=W|^b^4p=DM$-u~_nFFc@u#!f`|T0Y;yATB7hfP&T?D_QCWlXCciX4v>7n6IBu+oZ$0dN3e8B38wu9?@kAp$()z|%fY0T1GH&yX0u{5CIA~yOT8FW4M+x)>!SPlF$9sEio^fkZ#{fGu!0r>MW3jAHMXK+2AzTI>WTp^^gZrCC|KlPRw*il4l>UrV!!jGU^KUL6ZOyn13Mfi(<burl8qrmJoYA%t$-snbE<txVTA|dNJ!|#pi0d)EOvv!kabF7YwF~cbYjBZL|DOGovthCO~4`uX6u{<KGh=H7AV@pg*Jx`bi~Cj0!2xrz*rJe_N`N%FD=8hMrmO&QF33%=;A^;(zzsTu+~g~Bj9`KQK-GV+OjiW45CJc>JC{X8;ONF=HVqsMJevOLJ;;g0yGr+zzI5RT<frg<MNkyiuAMsw%Uq@!2e?W2_*iZbK)(P@jk9-TCj%AiI;Sc=!Iwv69f-!N1K2%U=ZSjQ778qfFck^5Cx`;r$NM98$8sYPveQwfMTAQT3RfV$fsEnSCq9rc_Iteh>j`%ic`&?H?=BvoXI7h9;A6hp$Xs^L>@v2BnDUK5>QWIu*Vr-@8dqNf}nA}6Wr5-{Bg&8cQyjJ8HWf4%NK*l<WYFoq0yVr7gO}{dgxi_JzN4(zuA)+s0?iejX%Jo5pO*OB;0HMwPV#Zb4SAr$|ZvU4-)8~KqMbdbX}S>ntUEE0ZO`W#=>^k%m`F$b1j=Z4$n~cZswEV|D*;R4)W(g(nvi11R)85P6^?9#TcQqZ}AKOD2Zl(&`|S9!LR`BIWXSKT7zAiU<s!Bk}%X&da12O9~=p0RyjQYLd%1Il<k3oK||||+SpBg8_$FU1Zb?i3~Iu+w$`|CTpdmF#h@~a;DZGm3r9mv4YXaL-0l>RY*8Rxg%NQlyFl)VX#aTF!ZhH2y~9(z1dJfvEMS9y0{_&chNgVJVivMTbn-~0Ktc!M4rXMOwAcvSp@6gy)zooL<|XWeIN(0TmU$W^qmugt9t(Km6c9vhq}CSL&G}Bl2|-{Z7_WVpr|+;4d)+}uZ=?d{Pb}gykcVl8e8IF8J0+%PqP*w<WX%HuG)}89(A^R1lrvAg{B|h76_Ly6(Buq*DPU}Fwnd=_dm2=R>C|#}vHXI0W!UTNNTsnoCgZF!M7D$%DKVv2h3*Yy9y?`*wUf|kP!g}8Q<~)+IL91t+4AU^Y_|AZuV;{*hS@ZCZ!zF~IEo5M+`8q1QVr*yo0|PdICKeHe9hy<paLNz`3L~s>orA!QKQi?c}v%gWpIU<?hh1j%)|Ad`I3&{6d9Ks_wh8i-trNUC3@$mDnO!RB6hcY&|(t}q$CUm4=9ob`kfbe?;jO^159%Y6Qec!AYy>M90bcYm@eI|9U-KS7W@dl70=C*5Jp5ScF?|(kXrr##<d@#KI|V<V&S;MC$S`iRBdV)==8t_aKHyeP9zde0ja$M9!Kd1b-V=HC<7yQxxFQ6p#n?*7SLl}z)TdC{5K?_&w~H~>OSMb^tv?+kQ29UfMT<30Xt>0#Ab~F<jvQEMp;Wp?Gg}3KsL4OipnH&vRgoW67b#8TvMEZlrLKpm<znmxH{x(M=}F7<ed39U3%sOK%MOUd<Cl+jDhy$7lR6k6|%Y7w24EVv^G+dQ21#n_YP7|k|}Elc@q1fE1|c6JmC@$OBVZZqvfG)SjYt^|4;`Gtic+NAyX~^NjN}OG_uh=x>ejt5`5Pi(cqM;b~~v|k$}XL%H%vKF9uOhiw`|fuk&Nz;&*<gg8tenhI|_}AhP-$vj9kIi~v@}bqk4>9Tbsib(Vx-@xrEKCIDh^p<RJ_5>IO$QA(B30}eKrBLQ$zfvd*&>K%3JhT?*(ouCSnfe@5-nezJ~@E7u_q}CQVfT)nXQScFAg=z3EcwEWLrr@|emLLz0j&Lv62KLE0gvQ+B;lIp!Q4ES|qWBI>loHL{F8TC(c~b4SV$eU7^!Qt6_69S?5I|@(OP!c4C^HOqaJP-8K?SB_67!s%m6m7yQZn6M8Rs432-ix^yzMRvq^$^lObm9OA|{ZfunflKE@3^-;4hE>9VZ=X^eE9T;t2VYk&q}?I+0K&rWDYB+=I>A#6zF4rVmpimIYLH!3{W%^kNX4tpy4L_;KSV@k2QgNGIpJPbnW{0`6dIIUg=GnI3hxM1;GIV31DB64CgQDWcSX6hsG@#)C?qgws^RNHu2$*ym%$JpmIBaH#FvMRGK(6aj|mFM*UE@F>i5_WMB~x8d`qm&XX}9<rkZOKn-K6d|4yAr8LZ5q$D%IA&yh$GAiUISF%hSgbSPF^{mg2Z*uh1P@L+l-(e&<0T`RO>C16){N=Sbdk9&>%NJYYrIFm1?Dainvat#IKCK!<3~jZC79dwBE#3rA^>OTe((uy1KXyq$p_6r*GkA^j5GQK6M?jBGcE0sk))nQ#DFmgoqa;^jXl)UD|G!j3HszDi|!RwB$S*JQta^*G4x66OBxUnd6dn}sfoN8r00_}mGEH7ba7#O2Fgg~XI{89M1%=|Ai3;INeC#)te1nXGna&{iteH;sI6}B5tfJ`J_!2o0lXf%c^9(|jkv&CxUO4o^E@IlVQ|gInWg{=nH@8m7B@s3@QM`PkEbETP9ZD=6SC@V$wLuMDNy5_CL6dms%g2=Opm%bI8DtA7%%xi(t;|M`d;#3M^9bjLUq{N_yG7c(V!3Vf%7C&4+|={3<m=W0n@52#}`639zTq;v|wDafe0M~N_a5{x_O+%W<<;U7z%!Hk7Nu%sIg^i4p2k!K5$%lf~Yh@?n~E5O#wo6AB8Dm0oXHj4aMoWT7#S1OsHJ0(gBiTfK~8mNN(KR9V$>(572*wpQzKzwAsT>a~`J5<h2ei>h)=ic8iJ!j?01*MZ!BBlMFL{n>52*=8!HV+DJv<g|PgRNqyKtC=&?br|ndfh+kp(k~!c}Xjlx<ap~(4x+fs;|3+X%dd2}43M^fRvUd;*lLDf{q>V(Ip~{K?5OnFD2kDPM@KKbIUM&gXnJ{w@q3kUoao_=WEBRjeO@s#^Rv=jZa0_9}m#slSV~I$xcS;DQWdX#aC4w#8j>f4<v3$Y97p@IvTQXbNVY2d#6^$nE%J9WiDz0@QOvZu9E)VX-pwbObq9QY_!vACHP<!tAid|TT1}fs|o;L2rt#cXxxu1!c4jDX@>T$ci1*Y(+XQ;o0yNzqu91w$g{FD&q2})_|t$Do@kqDqYgtT^^!t`Synk0Hgfvph%gZ39g3YHX^M5<R;ffLjd9pXu|OPz<53@#x~DtQUx%K`;#jbTkb44whxt0pA9Q*pWk#GRjt27kUX89o5e3)BFYq9F<3_#?4e%@#8e)UFE+cclcFljZ<n>o!F1kO24Ygmif^h$?Twq1ny()@vB-I<lp?jdF?zNQ}0a!LL~&8WG*aqdBF$L?l?5sMuS0^gyx5)X2ufQ#EU=h$mO-?YfA79<ynXGLYRH`LXJn=-!rhA`zw8U>I>@n@I_zLlUT7hwbFWAgBxAm7o<wU+9Ji(oT|7sg{nVN-@p?&pijgVQwh1hO5%!<1WI93csD4%_U;rxM-s$s8Ea6H5nwvID1Bbxdt6i%mCUVoq8%xq21^X%r6<MYr|liCEjY59GMiTP#QwK7*xn-149dDr<gJV;%cr>!ExF$rsy&zY^@2Mp=t)4D6i!uV=}`%F@PoFhMr>Kb)XDRc7Rc&<fNZBv4SxJerqagTBbp^iX}>Y(wP!H%qCHwRcq#^07Odng=mlCX;3ljki;!xcZ(FKV{K(^-hJ#T-Xf<Or-1m{BrhrAy#$1dBtco;mrVn4{EP$>Mbo*i^dW;86$+rbffvB0d8YKuki}^?XmQ2zfsflMH90f%q#=b%>AnjENpb-Lw@dOgh;lhFzLCGk8^=On4m_Y~^==_|0ji}i{Wy0J5Iv|3z+Ddl66YFn0rfU60XgrpI2JfKoXG$aD#Qdpy{<z_-h)J4^{u?QgQ3~oVh{QB5qU>(b%BUZ3EioTlD_b2kTsK)Emb^eNl1k7n1xrnv>-dA#&xDv8C^r!pxJ+QoUcP#Ok4>zd&6g{D=RG`C<$~mbdhKzr~A%RT43&R>p*naaQi8v5fg!g{a$2KE0CS=14t)60J~6ZVbebkqL-sl-*l5)@eqo_bd(7p#RXSDL>enW)=s@UB2v3c)QWPLc8S<28t5b}maYNPIHSRk>r?A6fba8Q9D_a~TSCg*e7xLe_eS=RSZ$AH^x4{!;sV$?qz)&M98~ABtWq&K&eI^epaH!YblsSa9B3%&FQX6Yc7n6{9EiJM{UBkm3ahTHDVKy~4?^c9T3#>+hbbRNK-MT^LDf+wm#u9937}yo^B|J(gj>nf3v$}iE=Q&!l2}2Y!q$A9Qc$xplHl9Zn1X4U++jKFKnAWY^bWFYYuFEBHhqCSAHeLEkElER7Bp(S#47+Qv_HpAiaJ3v^$ejIR}WI4KxRi;^WlspPZS%H!h_o|C1!N05>(@qkmqEWW^BHHGJ?rmV$(x0{c2P(tjxA{*kMT;w_p@>Qdql2D-#DLssk>{GJ<z%P_bbWTtecz1;`awNjYv{oN<U?N|y?CN;5tmOu@r0;#zhS8B?-y16q$|M|La9P+UU7CQ`Lrt}-$mGM)r18?@GWHHyY6ZmdUIdWF3~<wgA}4Q_frwGF-#1PR<uXSdzAk<Qe`4BksfvIMbH>9Su83TCFZ0c%gZw%~#d6F{2%drQgt0|#Q~M%5d*BrNghcVNrd@YtJ65(tjdAhJ)MRj?jUWAG&@B-E{qb?Ge14`5s;ADm8o8u6h>Ugn9HHO7mOWI@tPLbF<wJB8R&LQssWth*#v{ZIQAztoM2&Nm~I!G$>>n*rcPsude--jw^Mc7P2J`T&@BlCLBDS=5_3eD-2gf&B_E9FuJ;oI=M*^JAYV5mX`6G0unk<)<Y9hd*I0J-iywcqSc^oYrB*RGHhzLr)Vv+u*nZyz=o9(pDH9;W)S9)}aaq!imXt?8AyhM2>@k55Ok9QE3j)o^bkjHOdc$>%f`uMt!*Kp(eYhr;wo4h$Vo=!qb4TgoqXmb{`Hac66qBz|Q2msm)Y7ltw^Vld#|@aDxK#x{{E6QC`r2t<j8?0BJY(SFK?`96HP>bSL#CBap97z@I*T{do-a1yDQFZ8*(F7j78<pqTe{>@+Yv14wM@9A)MQ5~h2{`P1nt!htP|iYmM14a18yO;o%t>QR~-@j_@sr?5*tUr*RKAZyiJqgiK&nHQR8ED33CbgbgDW-3|mE=7B2#QTd;Mzgtn1;tZyrGhK4oKSE(+9e}2AZ``IeIhPPX_I>##N$CoGH=oYc#%uUvMZ0>_MEC*Y&e-y&<r2xI;zS4&rV_|-$+>C;i6wky>C7Xs0<Q8p+Tl1P-+AiyE;MDUyUkm@QJ0BcXI$TfN8b;pe|zz2{e`5lOrz?5^PweVEasGb$xYfSncjzS7P5s*a<eziTj&7p$|3?8mY^h0@36IC5oUxxB{>M)=yJi?Z*R<7smGNLdFZ}b`9hY(I_uQvBf$bD5i40kn|A-4o`M)oTrrrytbmU+ZKR-pd68cg=-%vdy_t+vGx+NE&8rNmc+1Iwh?L~Y>!JqQqCx5+<sUSfq<F#MgQ9G8xa$svLZ+cA>bHvPB9TI1S7A;xD_I>&IV@h#3%t~1#m(Ct}UPO*>qh=;6vSM5Sr^krs@)qeisML-1a$jv|tKIt!Z>IE+KTBpod&dBF*l4-Q@Oe&1FNbYo>|lej5+NgwVVlCcWdJpfiBF9pI^6jG_`b(A6urR%?mqgb+d7qwdEMwjpBuL>(Z838t&Vs?Vce62e6mjynyDmxwYzW{emUmxgcA&jw@*w1&Q1VwyD$pk35GVa=L@uanHwobH={#!Tk)j0N*LJ6jZY^QGl^el?0F?=)6|{k1GGI}R%@2wYTW)(xsxd`sFrWi#MW=sRi5uX==@9jF`B7%zv4OEE*d6g26g#OZVrd@-@^c~X;;v+ko#^7DqJSWNB@hG`oB-Q+}M??Bw!kRx_EXm$j{z>87v%0W`;@ackk2fOFdkixt}#U{ISku2^&*H)QnG`p_*L8&IFKw!GtIVl*h%ObBykp4MP+@K5uX@7|bJ`Yv~_j>6!n!HjR<NS=iZ;6PbLR{EYQdT##6D^$1@oEf*8ju1Z&bTil`5ujdIVf&jlPSP+-{b;gf*FybqXG4jjtCSrGi^f@_18<8(tkMVc@7gx4KnHKjs)n%{aGk^#c;k=QVl?+s5n1PT{G<0<HrCXr_lqI!`Qo!F}xau4~fa(XJl(rB%D@A8e7@8jYXchxe4NfI*D?Wp0D^D5fE-UnA9aARXm#qg`yrB79wLJC^E{!7Lx2D8%YB`c@mH}fiRwpH_vBLrkI)W+?NzY#@(%t;c^0`_hOXu_Vj1v6<*bq&VWXf>>Y9Ag<+2e9SDzx+10a4E<3QG<6$51e8KVctx`>1IQ0!B7=&$;s}~`=$<EhpeM4acogg6J?de9lnrLbto%?W_sqU9ux|#Zeh3tLiBgi<t804nQh*R8VbIF+8RtVtL$;{%z5B!uWPSXZcsEQv`uda|#I1%U_v>#s2qV7bBixQ<LM_#L|_bHrf?!L03a8~|s6Ji_6VkQFW>>t#0dEZ6ouAo@;9J~Du8r~&pwrOnfVpO3vk=`TZWLA+Eh`~X4vJ5u>*-}#B3>Jj-a|D1&3j5LF3UB~MlitKW-ZzrIo%1?@Lym{*5TWVDb8nZ7kT`-W9%enbQRgumzf;Au-who6HaoE6lhEA3En+Cb))9R*s*wNauozNdnwN|iF(e?dQg*3c26kYETk23g!I|E9Ti_FQvi=m3?dW$R+)K!@nE1#ELce>M;0oSCpu$(6idKEiq;5LX0nQcDWaPvkG{0|$0}Y9<GK_7);eaUXa0I-(97R(Zp3j3W(F+|7HiQ5ae`w>}NC@W*=PAyix^;LbGzYl$5)!CEhAd&5aS(hi4v|6pjiTaY&f#H>#PT2|om=n;DVEXoqQ)C!Wv5{4bUHwy)A4YS`0j)@xi^DIQF}F}7=ghX*bIp}HaaA_H7niOw{nFvgOCX*Z-j)k8oCpDw$w)<BX}dYjF)!h!Y@PLt>e)L<Qeg{Iq45Vy1i|Jm8nxIK<6U5{+M=s+|P2M&UcU;o=pG_V|_^-NABeFpkhi1;xQCOjwNNCbar|*^64G}eSodqoiUjs(ZT4>x$Plfm5i{05@1&bqv%)gp2s=*Qt}C)UrghENeN)XkS4<U`BsR52cp#^f-Wgd9!!6_>}yHW^z7<L`sBSD<^BtVG#b&3r3CCilUR>u123b?{2UOOvzV~`5&I$CJ|4t$&!IbR$`Ui<s!cI<D~YIw8NDRA3v2%QRMIlDxNrcAUHgI`Nl6Dqm+ii*4NR|z>rLZDN05nYKwk>`*sD<mzyw)%s8-jSgA>tbQ7{jm+CH<@M5F+)DZrC&hOOH%K4IIs<Ro(gps_B(PU`m-wGV*jM8L**=KwbHAxaL-`@pmeCJMMgl|=NxGz0Nev2E5&CtxvYwL{8Orgro4ych)^FOg2*p~tcJGMf64)JR%EZVM^<fL^<m>Y8?UOb<(x6hU|p5-(%We{n^U`S^tmX|q81nE4X3vb1oOXF%stP2WZuS$YDr41s={!RW!*ZWc8GAo(4Nvgt~RPco&~=&wcrvE^_cF5X#Ux?YfWTKJ{gCa5G$2|)(bN52z1!8)RcZ@Pkz37F&v6ku#In|(ZsS%@IfA_Y7ZvF|u8=3I~-c#3e^QE_wQy1DTqalyK!td2io80%3SZtoe~TiFPr{%TaQsSK43N!G1^2yPI~Jlyg1l#-5BgfKiwV!uSud5^Jt;uO*)<JoL!OGtgpN=}0}pwVKn8!|ZaIvo#tnVbS8Bnw-~%?FG+Y7~{8#VKWREzE#psuwpG2ICVU{2_PyVo;$|7hM$`oJ{ADh3Xu9C!XDPWKj=X!|6=WH3}*YXi(CFn0VaHseX$Y{l+cYAAFetW5W#A`Yo?la0fs|rw0HTWvwrEJ2Ws{by*rq%<QU4+BP~@F5`l$cetz?Azq9s#br9lj`AW@)@Hgx1G=Y-GjLu_d;yrK=7sAZi~{$AUxY1Zq)y9*Y*=!Jy6mcIR+zR*tli{iTnzw#eL&|yy1SAlbxJNd%Pk$hk1WoUK9QIvZ6W0L-lIDtkocet+Vd#Dm4uz#C3P9=hZ#*xnC~X!TS}OaAmA($A{B7BNp28|eGn3RC_{O<H<MPzQ~y-^3Ru<FVA?l`=Ps9oYUE6m4W~FmZW3Hg53k%P+wY>z(OC6?u9KMI?u|;p$*NbQNRSsQgmC>AM(1!>VM-6-SC6Ndv<Cu>)~C7dmSQ>=bY*V!lG9D^_32}AE<wo>7pS7T_@Dk12CWtYF|Ins!@v{HZHLvX-6nfI%@9Z#-_Poq0e+e(k&&mtIUrI31F(_48fCn5itgh~z5bMQvTIUzN6>N0sWS#2nxTcLAaLfxyke#X+hRJ+Nhx+-Vp6HQS!F6Mu^|`zWH=VL{Q)I!7#|)HAzr;XG&Q0Mi*}ytdMA(zq#d{`6%8aOI!z)I+{foJ#O^v^IVthQ&pClz5U@u$zstFmPzZ8oa&bYa+xAh%%z`|)UgXap>`$xFlS`Oa63OMS<Om~{g?TXG3Ve{!Rm0R3;IG{%XL8M}@9OrVDaOg}`#=}%DP=H`Fj%hoY80Rdmv-n;POB~poB<^59Nm7V394(~$eBf>ZCs{AdvkXqe-M-HkMs`STiI4O&BH9r==8t`z{Hb+dqr<4yWV5Y3$z9>`B9apU*ndCQ%Wi`V2vcjO8ixlz>*B<C{E_9QN`67a9F}|v#dN#9lS7pQjRAki5x(K0fBbspD2BcZV?A_S-f&ua6~|;PnqBnle(WJC&?EM{XFu(J>))LkEU5q2rSY&gifct+y?em6Rb@+-98yMVES~RXQ6X;eF8zY{%RD>AK(h1G;al=kTXywn9XQ?J|>v#R7YZqR;aashMRu)0k=a97KIOl#HW<O{mG`3$sj_^7P}^yI-zvIAcM8jK&z16xFx&&fW36o(|$8Cbz`FdX}+1|fHwlRLxNM~)hODT5le@Awe)>v^?wn_A@mcipJLMDWSCM6Ips4aw_T3e#fn=_(uJf9&fTQX%YqZ@e7pY}`E)7;>AO^p>p3AS>aYYQkDNbZZj<IfM&~x$0|{UX-Aqd=Inj0yGid|K#l0BCBSq7yND*pZ&q;+olvn>S3gncNJ3%1Cc~n|C9;RS~3gkh|1cCG&4`GQ3E6?WiAZE}gO$WCNOFTf|OHA|jROd+n@<%5Ta3FDT_n{4r%SH2p%WIV3L(&xZptNXsIf_nMrVS^CUY0!zl^BI*)~%;;(A<Bmn>H0yOTS|Xof=yCa>1gVakprTOHN-#i%AbnNbe>%!q(cOqAj$uq48vLtzL8qe=>Beo>_G_yd!lD*0}r7HnR@{kty}nD2^|z4<KoFvBBZOEOY^6RwtN;JiFh5Hu;ud4<g=M-NNQUN>bPbIP$)ny80V026s5)D4k1C%0PvGiAlW$;Z-O@^^!#=tRar01=TizfV|-%p}BcaGQ47NNdWQ4uLhYSmKL+P#dYnaxxpbxC7N`Go+;}ir$n8#S*3!T(?OiwOMj3O#t)y_y_c2n;FOZKbacKLv$xg(^JsdI(zTb+7!Z@ZMH)BS{fCEeA6a66Ku}0#q+uyh`tWk&yS^GlBR~Kjou)WrX&_x!QV5-RQ#}DaM8?UaQ9aJZqM(Tod64UK83&M4$hstZvG<bb+Xdi+Gq^Y{3&^g%v46;wWH>5ajm&A}iWavDJ1=H5^s5?+MRqt>&Q)<QT=D@gMim#3QV&Td(g%zb#|yaQI73~{v{4N_o&!}pIwaER5cRm3vjrIJY|A*1@pY{5bKd=C-FmE%+)WRt`anZQqNmdtlLSzS;a6dBpGP*EDi;=t&R)t+;jPPX{iXD3kSVkHWJ&_s-pj)oYw%hnZ4PuvNqD=`q9h_^nXLkV3C+rZEhV@D;(cu&S?>7?-}K4i1sNV_Sq@P@qwIQhAc8hDUC#p?4Tqf}@@@HiOAOn_MT0c!1|e1>-XA7r>eVPG2?K~G0wn8V0*``?n~uPB(~o2zk)(0v#U!T#DTcWm9-mB#x+6SC3riTI>t0vzh5$UivIAwteYgqaDP}ZFT7p(D!+t>fZfF3hXxW|#<rwZNOiNEm$(Om2c5@lei$R4j&Tyx9AT3ME*%^5(L*SGW{YNqg#H8G1eA?kgp)J-2G5JjROt%FTX-cs+qPg&6^jZ4}%f}cGNISZxO+wb0+K9ov7<5o1Vb)_pRJFmm!2>8c^=<Y`vm3pYSEHz`kYxrLI^+pUPMN51KvPzE@{7W}Se%yEM$hO74k-=*xZl$QqYci1oN=F7P>3raO6q~wVz37WURf6C@o+Wi!7|)V8j`RZwG|+?K@_!i%1IkP3dDf+Nvui;qCv?G_hNZ9s(|XDmzQU^_p`bb2v33<<x@@?z0g=RO~4$zG$?}`%hFqNLNdT<sc|1$XTz6kkzA?~?cd<*cT;<9G2M`?0cF);zC~jIDcqn9!alTz`wiF^pgdkr@I3)i(qV_#@6{;4oNiMEtI<kHgD^SD<ci>>*;2apm0H<jQmnZzsC$IPqdrIpUxa*up0|``_mGP=fPlCWG;nC|cz8@8Jum}Nj4tgV4h4h_XRvGq1dkVO`sm#>f`uGZ<51%P$MTC&klH63Z|1vQLPo>N&`1tAFCn35GYG0r>?AWG=St3{JqQ^M)&>B-gfva(QXW*P1e&7R1;cH+`Z1V009snpyORO#?k*GT{)`O}Ib1hlrcG@EtW|fHy+O6j33xH8xQ#b%_MnH;64SKY(d>Lz?L@c9O;T`XJqi?V=3D|i{~#tyn(q`dE-^b)VQY$Xp3&K=Ga^V%GR~GZ!<byu1O$n{`$R!<qs)beTMA-hC`|t7WUn)8y8jf?z?Mz3g#K#Os^+zm9#@Ek)P+g%lHnLzMhdKf4Z^w)js>RzMh)6?ULwNnJM`J3ZzYSzu*##z4FN~=PDU6n=+Zq!#9Aa|SpA?_NItIM)zfHm#Vkp9f_uS#(210NPOmJmza5Oncr}VD3Eg!twmQ}8z(S(SiURF9i4P%MfVHt3<vb~5X-x$;^h014g-ZKb(#uoIM8aw@Ff{k;q0a!zL|0v1F^6Wm-F^&qI@gouMoR7uT#-?;TE-8-8Y;0t_XeU|U4C3VW_=p<120CEXo4YHqU;VEhi#(*Jvw>AZ5BvNL!S!f%~keL+ySJ3QD!_SSbzct6sqeTNm6xdC-JNX|2|9ZqPz5lDPwkS2fY*Aa)s(;ss?CeLezZSN$o$yu?F{2U>^CK>@o^`=l*IGmCMdQV5;>KKh<$l;bm#7ySarlAu@gQ{wo<FUjz@cGbV28yQ0sLsv`D(OhyTFx6M`?od)D2>FG1x+q{%ClA=EVJaVh76S^-71u{u<oJr|<5Rt{rJ)ehQ7?7c+A8;vq7YWTmJwoftGl?m~n8T{~-%wwG{#<)-Xb<<1_DCV>@GOl3nqyb;mDzx5ipNRA%|$vn4f=X8?G!~GQ#yJ>K^cX{*jyuSr8K};LIzJuh@pFMSn)ym_<l?QNy%G6YA<g*N+CYxB3_VOr+C0{Ok^_I;>q7(m<;NzPVr!HCFZ{HA+*z3`6PFFu4I2~MUaN*AmpcjP&VTFk{1|DK&lI))2m~g7cn@A7Th9GJ2>*4cpFx5Hq+JfC^#Z$fYhD-TSh?g3rng3tM`03NZ?fnHK4J~W~_l&fHL3q@Nf{{i*y<-`zqGDXrZ~51JL2&ZKN3mxtTb{oRu{Na7It_E+*4Y!6mb^59sA+!_<B2bWCY1L>*(I_1*Ij=!>Be!H}%0y!UVxVQ?dG9#36EvfXDh*6B1)T2)BL+jx*r)?#lo>Dn9FR^QjdCKe=K&UF#MgiSL{0a-q1A1P<y9YR2`D2^D*q+d?6NqzxTmuL(<A;QC^2^Bi^<)~tgd^#5qk6)o?78_ns>R&HM%UmK2s^>EOBw$k=+VXkgwv@CyA9@A)sclS^Kc}b*juChCY0mLTmx(E;(Aki&tnP4DEhh99)tvsdwEJcp!O(;RkZ_j746Hb~h(_7>t3lj5dY)yIEA7*YOB!h_1-_@8QV!B}U@4rBqZ264wl(>ql$b8TV67=7wO?VA#G7QR%p#Fy!lKDQQ|YZ)-F1@=E0!_1kuil{am0Nki&ZLI)Ez*kHjvX=&2#B1`8*1Lk>7w>Kj(d50Z3%fuDk1So;GxHD+|4$aEj_Z0?7@La!U^~Qfg7K(^hI09?tbTQeQi0>rfooZ1j|Hz>e%6kyj$|acLb$F?h_V?%G%aVjD10bWqa$CfDh1l|c0g<@7HG@kAz5+tKKJ3rMLwEdtY4`j)VC137tLbke%nRjZRdM}Cmdz0B{hxh)AJO=zZs90pKA=1u6sVVEr=^$ZlhksdEgNI37q>VmkuNbm;9jvvY{-E>)4PS7q?+r1i8%vnQkngYUeOKGYb$Z|%Cmz2~YROhJdHvn{@e<%?l&T5!a!r!4A2|#E`nRb^J*!S(mnTS4sg3|ZYDI{(<*$sGfz~hx0EKd}x0lKAq!n)Mv`s?rm!xJd8kbcn+L4xVcuSP-r0nUW0k;1kC9D<-(F*xz_{eX@1Fi2sRmZ6dfJpsvi5V9Ez25I0D(hba~l!dDN2rvf<jNz`)V{N4|N}aL=DG6^SNN*aZl7)oYW=Tn-3G%iSVT@a*?PGLzc*;D}t5F4@4uLvah;GjV`d*R~F0G2EDMz<McW#Nxy+Ryn^pn}kwv-%W+J29{mDGSNwNr5Phz{tL)#^%B9NoS$se=?<)Y@R0w~#al$Oc*Kdey7H3`N_ZIh(NtR5Vam5NUJ}UX3d5sElz#Gp;S9&O6lD*=0cIfvwQ<p-?A`;3FaN_y;Yrx85QCz6;;lY{2FN??y3YmQn!AlrblnDb`3<uI{A;G!nw``?^i+B%EWkq=RFfoy=q%^(4ql39m-!gPVzNx&zot$mUQg>}q3N?iE$u_!W`I++-cmP0Vz!$itmd61LaQNPR6S397AYnT9bJlZT;UouxX-66Mr`<)xdpzHE3e-9_C`o0QUh(X?SM5l@=^bK@=h4bcVl4bQIz(f)@9!6eVdwg$jy0OEpm|E#CzZcjnGhH4+uy-dIWL2vh{dZ5k$y|7S}v-SIuk<#qmW@yM{;gqIM8Qne<$|nX4WQ5=d*eB19PTsoMT$-br+>r-c(%k8R)98^-b>)ju#g)Nuymc^STgVuEX1SU9d?;xQ%;O1i2n09gLZGfV&qDPnWtffb?{uY~QWp1qn^UBM_`t}YC!PkbW*sJ^4vS0Ywa>F#I0a6_0Df%BazA+xP!%M0&E|g7n%o`Hh6NH7^@~xYLum^U(0&^c8#(vMfnQfjZ!ss2M7AH9z)C>1L3M8%^+8NLZOGi#wJe5Y&KU-XGl9w5KqpQaR}XDk=xKp4I`K7qU@HtgvuW(g6bFA}CzLNcr9ksV<gw2lz)u6HO(FWk<1qS(p>{niwGn;v(}#p>hcn;=fv!&;?5JJx^tg%O%jr+hmuf|~tdO!6vT#_l?_<R;3di(6s6eInIc<T4FSx{`>ITFZS>*dJX`v*B`XER&_1?mUKIsVz*lVvw(ZZBB-Ta+lzaMg9UF@0+wBAyJzJo0zJM7J8fW+6knd{>odWlde%gcgfVFi>;O99Wsj7igX2NEvtV!V?K?RIrQvLYs!d`D|Q0NvYKsz|G)VVnvkr&PQNy=V^u`P{2f#Y`p?z0y!_pVdL{#!Sj!lzeGd3hT@P;-psa=Mvy+GY`AvG=k(Jq6O2qF}A4r>Q-;q6Xw>{f=(%E1P28}GSV0+sRx+>>qPH^C8YsSx{p3d@kWOb_{nia>FWM!6eMT^7BeT(_ioy-1aV^mhEqy6;WNy`QdYB07ptcSDQRSTkmlYNIu_y5jR%F9C}0Y@Tr33%!_4=@dp6TLm9aELvz*xt%z4VdmaKD+{Y#%*z=cUKqprUzPNKPYp2vWUfn<qqiszAC?aC63YL3eqU{CyE6NyhVZ~252(YJ>qD+R7!yr+$7M)Rd@|5sYWkD0l{oA7C#A_g*Cfc4ONR<1B{aqbAv0>fG|0!U5hEZRuRj-!IZRDk&U?bRsPez*w&2XFaW)SFP{1Qzfi{Lvz3mH>XX+iaW^EUw@O9;5_{*{Px3Tge4lZ*U<(L<Q*dX5cZWAR9Sp&7!9cii?~~7A`p<#-p5(p#6PgaVYyosS`#dy+&6fii!X>Q)cjL82yrX%w~<=lCv0soy<aT`pr;;Wp{{4_U)m019VO;{df=)GboY?TK5vT!ZwY;pVRDTLh<3iF@BYsNe3c&#}h7~LL4z8q!31`NHO2Xi7lr&AI&TvzQI-2;Q`{~77L`}pNG+i41}57;1o7~n-CXN=dOY}@sr(jhlF0^n8_NrYY6E>nj}GwB`DB!SzL9?NrUeB!DNCU3k{nqqdB}b#RSDk!yw{-!R`q7CA1O&)V9^^;TBk+Nr|YK)ULsD58F_rWaGuCOwJZ6_iT%LNg3T}ufltm;`wM&;)n-I3fguHjOJ8!N_9O5$-CWYfxG<YiZ8HZo>1M;s35FV2vt8bt3`n-xmQy^_L=XA?s^!V>N>WqZe&%Y^wnI<foL+>JLm|h!p^UT@wr_t9WKl`Pc04QGCmxwd(WU2v!Y)zCxsCDq|t{Ua|Q~QpgdN%wQ(myD7<Lf50K44H>ug-CGrlFmz_Bu!lAT4HMqxtSvU|53SvQc-p?yeBO{6h8)+9vBs&yp?kS+#8ea{gjcs-NtW<Z$>ci1hdW4dKIiF3;q65N>GFc!8$(%y4)jXM7OzJe@cUc{)>FCq)=?J(yBra$)bs%*XJ<*yO8Vt@7#r>9&cNq>DFjw2ZB>KXp+TM%RG&YI1Xt-#V0>I~C(3TeT5~LM*KPecT1Yy~z6Qu?8puw_f0OfxnC(wT|5If-vo9gMinhB8dyou<KkdP$v^;E#-bbh2)`<P@Cgkan&K0;=a&pVzRVxv&2(0GUkk|HG96vD#yVi?ctAjF%^F}r3%Js4SfIPp%V7xPflSH>ad3~JG=Z;?XuM?nF<xBCec1v|4*E$C{xajj@B35y(8yujE$8q~uLJ9A!MT}V%07dp=PB0+nuGVBDPd?|SLkpzvNUJV0L!S*xdwOX)-Hc+5#IyhTM^HK!Bf46`B=AV~HnG*jyq+?5_VbW33t59ef17VmjP;LU^ZQ29fHK@UaF_Sy9kZT728O)s8&&kmG7!(S05jWPBYEdT5+xcpk&IcqL@NA5Q@^Fllht?e}Pt$Hj;1NSd`4TFUYC*BBjUpHx>|M+{c-t?4xXis5*VB)L6nYchoLzzg4PTvUtKm-2q?wayalkdO6P!IQ^_|2phI0(hqYihDppB#tU^zk3z{kL2z*OB&k!WoJaz!`#5(!!=DB2h2m#naAasjDPOi{?AEmU$n9)%s2$~+*Z^D{>Hur~Qb$mVeO+RDbtV4M?*CWC;I4l`7SF6a)D!fldwz8J|2=U}d|`*cw8b8Qf`p3UTN{$2gE1aM2X^P90XV1k+K-Gs=zO=cy6+d9dvG&aE-^NAoMsVZoGfa|Q!iQV;Ju2H2g_$L~*-mCGa5QD`RBOfF%hjIgioh}K8cVL2oh<YnQcA4WEOY8x(zC*vHNmqNnPiCY7`(U8m?0)-+-f9#cb1mtoygtTdvKXh63Rj;+bfUXQm<cUkjJ)P=jz74Yaxcd<e_r>)iKTb*FJ_BeD;FH0pfO7%C#$F2Jgl9&4?ll|mks_MBTPQINEhdoN~C+;0zi>ekEzprU<f%-Dq4sTk~e!Ec?m$LAIhD;VQt5TdzC3`p}yno{AmYX_xx$3jVFki0QUhdNA&zjPV3gv@uHGno%7GY(cFS1N|){f1opDIqX`^y*rf7#VGMPlwu`x+Mqcyhbqc&inE=Zv0O1WYh_rVCpm6dQfCa<eq|<^4^I^O63F;17I_jArnw$ba{*m2NkH+<7eT?L!vXWXL2tD+By49Q$3Z*CpdNJ}6fR}`J<ZOOEPY~~}!8Y$Vz$vCZe|VugNppf!wYhss@dgY!+yH>#j+Z@uv~(lYMO%+O|3LRT^J(fT@R)WAx#4)A;Shz1&jXp$aHB+`1x6BD7p~bu3OS5ge|P>n<fLo%pxpzQA4dpaJ%KyRs8v_4WoI1=l?n6lpvIYWJY-|xJU`gchE5S-#zob%schr0(9i7EK<3oxVnqCrp22mb$LUiyNq>QNhX)He?%SIj^eIFy`+b%lxLeLxMX{8)Q#E`{2+^mAInrmbliu>NW{+PotjmTw(x>CzorC59$}fg7qvl4DjgEN~e%6P%69XjnPPXpt<Rd8ZJ$d{@*ik|9)AFmk3(_Qbq$&W+Gc9iFoO`6+I+;I{5Yw7Fjeru^1;{FiX35SVp+l53&7K192hQ1(cYso<n|$MO_5{`v^vLg*9&cseTQGK5WRrBOOk8j{=TC?pl-^)|bN<Or1*|Ugdf5jbW?`CEVPo-f5w)eKF~DWwTrUQ$)h}etlF4t&+r^L02QAa+ntsM3JEzZ~CuZZAJaGA=Pf#e+ug-Wxtdw*5Bx`Ztq+;)RPrp0X0JW$z6w`;RNai7~BpH_%1J}XNj_T1pzJ)jEe$Zy)0J{Hn`t=y+bD3RGJrfGp*dS#0+G8rhZEP@6pOdFj4YCK~&v(+<bM2#hTxn&fHs<UMN)+K{H66^0fphW&dQGAQ6l}bid~zpL9k_RsCj`82Ej<g!#G*aulz33M-1?q=g?gmf7Mlv~lXY{aHK}HMKi22am4ugzR&ou@OrgOtF9r2Q2>CQ_FFpOo$&`|fc&l@$z#<RvZu&AT(!RIZIv83smk~xfF|}EP3_MxapaT-ARMF1qH=tgi3H34et-Q$RQcVgIcx5Lsc@DEdK|Bg#;A!BTzZpGT7EhJq<IyJy4U2GoKYu*JeRIpi_#2Y|ZUaD2KGHUuL6=2pjG+0S^B?Yb;@&IezCR}tnIT^j`^11`m<AlVq-%VBF>X(v5@hmQwa<LjHGSG15MFyTdzvA&)6%l<(p-Sh=y6E1*|BFY+mzo8?pnG$g9v(+;yY66CGFK8r11D?_v!mY4|!6q0o5J3%1;C5gDdNCA*~?DfAiqNeS@3yZuMsu!!><SKZ8V(eHh4C%Q101Zd*&^2LRY~#hTMM7>lHqet1uN`fefDoug>>hw^|~8)Sb;-s`7<59SUuY($&L37ytXHDLizuXoc2ini}A<?8Sul-22^*~wVbM`efEW;plx;zG?yA{`^IomV$F*9~<Slssr20*}D{$=yMP{ygxSJ}((PyWKnPT>eCapf$&i-_GC9^q#+~QUQsfI$d!T#;hzein;js!v~R<07N?}$)sJ^g8-0`<5hzuA|`)C`Yv;L%?I#e;Dfyk+R`^L$8c^yqFBHLJntN?LHKMRpbCC~;x4~iGW@~GmOXzn;ZPrqMLXvYBTIKaFqL~aO&%bh$rUfgLDXOf7xpE6!19ZM5B>^Nii9sqyXxV<P8GZ@z$0(xUuVyeIe((xgp&=9m*vX&taS;1qUQF(jlGmV<zr5zK+bh1+&wqofacDm1a#1p0CJeICV{-8j-JLn%)gtmkJRJ|Mj(lFIK1)R-x}Qv`qNK>)EX|>VqBlAFXRr`UXFJ=aJT|+)y?ZGRZPBwx46lp-S+YzA!hg*c^#~tjI)zn-<GUR-w%B7ck!`_XTxQFo_~_Oo#uD%EoRtz@)iIvFDOuO<g|qZkGb3}T4cZn!}|<EZq6UL5WR#j^~iHf^+VdA<cWKPLb=dTmKp_^0NA6aaSyBS;1jsN)XT&2H!}}Q<~P$%E~eermqt@sVfeU{S%OVWyQU8p26RZI<?EcjwgQ!4Frke4^XWPOvFNxEgMxCd73j(-L2!r57Xu&WPd2$gDGeOwI{mo%VAs;S`J<Hf?Ip8ENvfVRmxKcCI+q`Ac|z8tmW~;JEJl0Nu~85EeVZUCSmLfRxGYp>h5FqAI2>5Pi(!wmhoi=nSA9I(U6u{+{r&1=h3&gbH4M!90m$2+MRQ+GplDvF#nbZ(9DS^-Ta0lAlW{a?CD!aoI3%~NnX$aA`)T5-YD%5_V&KF4p&mQjX*<O$)<7*iojbyKH-GK#nm+A)xptk-qhvS2M3PJf@I;|U099P3J&sDmK?C&Hb-322m~jenPnxSjbf3WlaT<tK&=dH1;G5Og8q&PETt^=s^VHnQ?A_#H{n{_i;kO0U;ABHRz*4N)yX=)H2T#{D+AzB)Yo%OD&ffeG5X0bI0+-J~##%tuP(>a3tX+LK?%~ynrklt9ttS(nhEfhY;B0<3c_-Y?$@8pJ5GorT;$kQu&dJmG4`T}y;GX;e>6|Nd?qfMyZm(-v927XS2m8@XL%@6eYT$#_0658@x)x8$VrU{EM4?IiyXlk5*-myagEnn}06h>-YqLm_o`k@@KrO<zzt+AOM-3N&)^~jfg_V%ro4eL6sG}TB8Dcu}RFCB6agPD22q4@w5GN`C21JT&h;+Z3Jk{O%*b@VTa)CB3y5Ue~k-6xeycTj6!$5UEClA^bLUFhh?qSV7J2SyjL99T7@J<6`oi^%!?fO0qyrvI6g2<<Bfiu2FUr8!9YI=8cNCNJY14h2lmtu8WNHDGTHGNXHxP0vR1Z#eBE5$?r<Dg5)I`{?tEI8EY-oY`4S|-m%Wk#sI7&yI^2!4Qk9T&g8^;V!XaF6J>v&TN#w^rMRYMI~^pnquw$UXhuOzHx5MV5f453Y;uuv*s!7{;Gs+6stP1}zZ?&4Eff9Pv*HR^q2|%t#UvQqQ)6VjPc?ACMUu!ufXY&=u^ZcQ%C0p$%|39MD;Z-j&htT@c-ycNQ!OF=Y3NxNm!pG;DE*ikQp{W8O1FjMd$#`ip_gsL<1yzG0x5kNg-JBvu{ZynYv5BkYjRxR{}vKzZEK(ZTD0No;~yz&g$OZ6ITF`n>pNBi3~9VeVgmIS+~p0^`solB{#(BQ$PQ<UfyNN*TABY9w%GOR4SIo41(tc=uA_2k!6`FmVzUAW!a){!qW>rO$WwI|f~RjIeC{1||%41Tb;;G2dBq3VIca<6@G(40u%sN%y6?=HhAGoI8jBH}D1N@@C~_-}3y*MtV2*dfr#;RY-kI!a`3qiYxZiIiTD(tjDsZ4zns05Q9?VbLNIpQwiMG#0l#psJ(#fy&tvbjkcg%oxsVvc|!uOWqUVo5~k->?Jg9F{v<+6=C7Ld|HC*rh++2=1rD2v5FoW-2Xz8<&_)g=gSkng5aK2e=xp87s5$dMms{f1iQwzQ%;`sW@Z0_UStsUrUQ{9qN1iAfpsE!6rM2h>N0ab6-f~-Ko##Z?MP?*;?zQSvZ)fmoT2#(@hD#P(rp(?c&*SFoV|Td~6gfwIb6)7N0}3O*SCLe;oi`QiBCx56lz`HvO@tw{CkSl1_*9oo4<yPO4X*4RNXBEVxV`W|@kxo&6XUN)NxvAygrMft`KJy%&}Qv%<)f<w@s-}rJ%$-*UvnqlMoH&gC%fZ8gr(eZKcV8HS6f~)haa^=v!X+KFVmTFr2^#8FzsTK$+EjAw|B_K1HTx@guuTQz<Qt^((xDqCvO^O75clm14B8FyQc1LM?{Bg0<JhzCc{m-n{MPs@9WGoTF?`$YG>ldFmx$b6RQ!&p)Ds7wdPglyPn25w?7G<uosMrf0(<mKZWKZ?@l}nFK74#_Z@m(=z(~mh{Nv#Pg5G60ys#<x~#3((K3%c&~`p0gA}2thDM)kM$d5C4q8|k`uAeooIE|AWeOL(pSP1w=H?%*#dnkE5k6I<iPYse-^?sR_I7JXYO$#aC$qH_By8XmV*+a39;`G1?d=N@Ce!)04sEg7(5cgWHEu;OiM)r&za4hFB2S&bA+vus@dEkVo;ZA?Fc8z%iF$AfJlD+Wo&W%fjOVflxQ*^GU17xN^P&PCAigY<kq^{w6`C)CJmNg_0zQqKs}8IU{)le!o;Rybf43vNdz7N9`S;ZIFOY8N%EM{f#`&I3`aDN-`^*M>>d@5VB~xj%(<Lwin6U=%J0@iCUChAfBaze><K{gryi05|4XhsH_1z&{@-6T5BhxvLI!QN}8Ql$*V10;t>e=v9Osuw4J9HT2tpVd6_qLQ{<HxWE1PBO7j^ndJ25J8D_+s2tkK{JZ8Db)m$Ej0k4~Ktz_oTvBk<YmoD!$x?_TFZp&$aF@$Of`S5OSS$uuI^k2*{+o^kXz>COm)+?4ciB$!v1v<JGvi?i7icOB%pjc(d*lT7$lQcVys%4d*q5&e~@tx;w0q6#JeMC!@6@7Ua?~EG|8j0m$B0#4tio=Vq?NCU~EZ;%EUicAd)f=W%oLWQSdOWv>1=%Doy>Ro;bJlUMV+s?<H8*knmlC`w&%UDgb-VZ2$m7PkcnW|Nk26bNh!l2UfaFaeoV12%!)B)9IukA5|7y-WrAU9gzMncqwv_&#B^cP~>)raa#(Pd|5Ryqwi1#Mrj!4jp+kH^qo6NXE+5?s`=k8pd_zt4Yv76GvayFxw?#9YG29V%(fOs0#{hs=KDVIrH=z2%Q?Bx0BDQo^!9zD6`jYRY|C$l(q25yz^u^vS}uLU-gaw4j~)1BONAyrsS4e2bKYRD`5kQ5gK%JPs8SQB{af?#stl8ECPB>w&d1|Z>LT>z4N-dr!u&2`iBGGTuEh3og%mcXe^-7rQYdIK-6@YBXF8wgTfrlE;8rXVGiok)fO`_JP(^QAG8E<N6BjacFjDPTp8;7cQdb(xO2{Y@PPCN2nT~T(?4)c9a!bTO{ElT8e-ruC?xj9+*1#;*9XmtO@s@hxU5Am+{pl|r(x6lb$C5NGRC()n!m(z+H6(%cIpsM?3<e5LybIzl<3GM?QNT3#e~~rl*!r7nUiuQ!a!$@GblNXyoa%p(7dru(owXJG_97G{r$k{WmED89G{1%K0c~;(z2VRzn%Lq0pFTCq#_u9H`Felp4(QD6_|IMxVUA@bgrbcHhOeV2n!DeIZ6j9Tn`~t0x%`ItcssJ4V>W(A6N)K-68LB@(@fw&fDHamMA6OXU3L)2l&8a1jaaweP-Z{1!qlsJ+GSzju?XxT;6@O6(E3vklr=G(OLLXaIzC6D+MN@=W%o4p<cj%r60%}Ekk4jM=ZTtcyr_RbMA%CUD9Gi;fQmbPZNQeMtH}S;cM!CItl&Nx)b<N>zjYheFAuzKzKk8P1wS?o~Lnh>QqPT?h!YO<1wma5`DTTzt4dHxVYARkj0_1DK%|G1LApTC^sY6ho?W*y2DOH6)p9P8?!|ZUfPWD8cfdVT<=&-M3Es<e=%-O-d#E}#;Q)v8)p?>M%?eG$pfO^x0Ns$g~;R(GVwU5@Al*eJh#P$aaofegcNfC!I-H;*zt`6A}E@6hc$Vwi~D32x9502ZcV<yiAdp^Ue`COhAtr4U;J+JkhARTTBz6S0O8vmm>=QeZH_?^L=qqxjJGBastglLdh8zd<hdtQbYTMe*tENOVX+aqv5W7=tz$zAcme|CiQnAT;=n6(Nbe?JFMwe*8B*PVb;A7<QjF6k_Y5|daja^mPjsolv2^jtR+DZU0UDYrnHim7NnPb77l5BUja%yu$5NfmcHSU{=<Nx_zr2q{LMeJimV8%%XqwEHa;uuco;)q$^jeby?^8{hEhxJRzJ?7E4fdKaN)ok`Og_SpQ3(H&cY85zPJY%mLOEH+c%x~myEMf2Xf0@Q@)=yBD8Pj>o;+ocaqN}HeWoX0aqr(b`A|&2*&T@9x+Z2L&&|ML;=b&HU1L!XtA?j>>v4s9o0e@5qW8Ctt4r<ZoBL++#RWop@^CW3_2%xb5O%a#JjUc<OXN%$oHr*A#1KP=TEzQ$=0Ky>0euqI8`1YP%mP-k4cYNx+?sqb7!JnV@&??>7C9)H-mSc$Z~MAJ1rD5!V-qIBW}1`NiKN}nAm9jqa87;_g2MR{sRLS*r|l{99u1rylh*{vdjPwG1@3v=o_sZi{vaj(cIA0nA(MGO`FMrj$tW4yF{=^BXsd#O`ZrjyP-2Aa=F~}L7)c=qaz_e)m3R0i9jfxWOUJ@}4TE^WNqHJKr*2ZGM0@8+#r!b!;o7tE4jW1!t)0OoSRm`{qJHmYqAq7mgWDOxi3Y+ctjW^=xw*Uwk9kf$^Mt~uP1!#>{{J*~21ymeFbv%@!%=>xN$WPb|HbfAoc!OWFf3R>K$PIxvMtG-kLIiN6}N|T!D8bAY+4X2*I0wq5i&Y<f9^Zw*Bks%v_z4(B#faJ>rb9k_>>B<%*q9R#~uVbw{=^B#b+)AM+)?0GL>w~l7>R>n~XR=oQpKREeWwSvR<hUs8TAo>BchyFm4UwS@J-@xk;`TDR=NpOYs7xml;pK$gGWqyeu-K?xSG{#Z{Xuh}V7WwNTUHyf7aNCfGERz`l_O<0}))(8M~mDWP@;!F15jl(KR~+riWEGv?wvnavE<6j@kX+<4wPhHT~D<cu0y0IR<-U`kVlf*U&1uk$W#%JOh7YF3Y)Nk)o&?y-Z{j_V->sw;Rz$~yL<yv?I#MV~k;mG9W~qgQpDVlV=G*}~#2+xJUpnXtb=cS#9AF@wFBiaHtNRWTjYxuAlZ&Nk}->HNe_rQuNmiMtzil_mqcmf-!{AwbiVtI`0Uv8FMEZyUS+lAj)%rv7DHZ*Q}{ut9882~fc(vG$2f1WRQ^wdq_?!K^JBLzpN?UPwu1i7?8&V)xsb2E}x!<}vLuYo6Kv4sqQv-PN2YiB-Qd6!E&*hFfVs*UB9Z8%167lfR)cGquTF>EYZBIvVfw5&A7k;u<&@CWskE9hH#rn%3G*XcuNiT|_!YB8*5l`wzT?J?&3{U2%pKK-1<xC;68Q8>Hxb><%4Gtp7Z9Nb&KnO4sv5tyCeBHT-ozKTuX=)G2r1C{JrFn3_Kaa4D|SRb#75$Ymh7B5jNf^L*}({Je~A|NhMPcUSn3*?_;z6=O{5=r&>us^?5_ma_@SyznmeU)0!8s!JY$_m9tfe0+cX087)7lK', '4.txt': 'c$|El-L53bb)Bz^0N;U^Ll|Am_>a&XHIxk*k}Q~neCI}E)FK<=^o*yQvgN1u`c_0_<>*wIfUW6sPE_U2h}getul?5_zJB+o_ka7l+b@6d`w!p$_~(!FzrX$A`|m$~{`y7#`G-G!`tDypfBf?0{Wsrz|LfavKjJTc@z)=|e0=}kKYsc2<?E+^`}FnS-+%MNpT2*`@53*D@w*>>{Pg*c{Q0}|>s!8$R#O_e$6jlVfAgP0Ilt@=|Ni60+vj&*^?!c(i~siH*H7Pl`1AX3fB5|6>xa)@Z{g0L{@tg4`S|tI_aEQ?{>P8sz5n}ff44um%#Zr7KfL<QZ$E$d^S|pE-~as&|NQaxUp{{R_|u0!fBL^4|M^$ngx`c;|MD0A{o}v9|Mugb|NQ>9KYsk_>tBA-cYOT%RS$m)zxw9yKK%5@k6*s(&;I(;hcCbSyH9`o)7RUfNB`kpzW(nI`rCi{_|un9KYYIZk8l6*?hn8J)xXvIae3Nr{_;g{V;tpP#)!GKqsFe+`NLO<`1?P7{P=afzz?6l`{?)gn=k)MudCdTkxMM~sJRBM`1_9^zUa~Z`T6@k?&AnuUlU@8Bc}O6K74&2%eU{p`R#AM`#$P1?y077^jyjjvi9Iz$iMu>-+cM{{;&V?=YM_w^Jjh+kDSLnw%C&XxK%y$$M5xmKKKoN^X;#1vE54wq15a4=sSP>r|*B9AEMvpd#o|Vo^#^w{{GYF>lVe5YRo;QP`VypzyI)mezt}8qqR^=?WK(nI5u_e#QT5z;nU}@?;n4v#eewe!;kNO_wn-=znk;N-~RB^PanVg`u}-%B^^huM-EwEGI)0%fAhaSfBK(szPtaamsj<Ma|$Cj{wy@TJUvNx_vJ62-~X*n@h5%!KmG9C`)~e<Pizvn93#Y}r$6+G`%nJ!e}4asUgCGUR=>XWdM`&zIu*T!nA5q-+`9LF(CgH>`2NRRh`%}VZ{iW_uR{5&qyH60=yxA}`d5A9pTB<i<41q%J6+1sLd*KV&USe(wG03Frysul@ZAre|M@T9&13faPk+>%@FJc4k3W6nhu4q#cK-2tQNRDh1L=c)`{R!~@x1j9KmGLS<4^zbUH{?3U;KrC`}yB=DRi&(@8ACAkNU&M7JT#NUHJ7)$KZ#bSoz|U{^J|n|Ni}7K1&1Xoo})J=<EDd7w|mh+Pw3>KYv=E_{(4D?HotxBkNGrp>thI4BD#La!Yyi*7c7tsxE6i%Aq5h%G2|A`-czz;Q#mZm)iT0LO*=lLpcAK{}-Y@SQo0-aj&JVXR5u`&O`Cb)SDlwj2_bw{77-6mb3JLZXt@aj-sEn(Df~KeW^Ytl&mwPAL-$Boo)rCQ*?=kE?*7l)zUU1GN;s%u0XoQbgwa%P(zDpl*+*xz1M?(w3KSETbh(h=u!@Pp=n9edKp#Hl#o&?QZ-1EOU*|ftszezCF`^FLDE7^Pt{8gdgo0SX0&k?Xz0_UG(*bzmRm{jt{)z~NTZ}4UM*!abVIe1x?Q2)V!WqTn%-WmJ;r2#QkJSbxXnk<MlXRvD=4dQZ$}GDqOND)&}996&|`B>^heUodb80I>GtZs^<*Jdwy5#BUGGcZa`cNtRnMizIHVR@9XGv}yPiQGmP6w8yja?XUTg|HRF$^4McG6hup-5&%?;k7-lKlc0WU%4a!Yg^sgKacJ=f8fM0)HJB~gk)T3OGGM0#6sa1bm}Yz4d1vOYBGv1)2u5uK5oukz?U%B0j$^@);a=<G-+oj)nSFrF`LQuSl#G;>4_9XNef=P{$CE4f=5HdtD@w=#~Fr!-;c$EC|=+-sh)=p7wgo|q(m?ZL-Q3-k<9tccP^j!{QWf7i8-**2@`MCvu?pKXp_ebnpYiUb`SogXPw?3S<JtCg+m0_k)|T|XJ%8>>Q!t?5AOF%U+7B8_0{AEhd`Fjae5C82J6FH0OfuH}&>()HJ^oJvVsCLQIi(YcnA&`+XXJ9i-K%;~u3P3hiWH%DJNq&2!urZ)Gcvv=1O=e+29Yksv<>Z$7wa_(|gLcWF0gR!4E=@MBWiCU^2_RvcYr$D;lx|y}h*Nf}+1igraKH!%2rWfqvB@G~RpRC__=>sdN6HvJVDQ4ZcQFKo9;8%h0q4ZAl*t)v6T9X7y5~u4Bmhft+RFW>6u75i);kVFq>1s!D>BRN9MS9YVk_!!8u_a1LZjg?T)Y+0pzxAsdElYKmKKDeowKgeR{!vFuC}ZufK2U#(GR^YIbl`Jo7kLI34E>_xR&P4J_mO4ft8|>!zFOJz-~&HpghRIE7OGB5mPM2IcgVX&9DRnaO0qYod%eZszDfm3ei^?dkkrrMhWFkZ?{Zj%qRlUICHxQh>#{3>0}y%3It;q)Nta$K`g-DY<fANWNRq0ho8E+8EEk>aI9@K5dXV|i&-ta)V$gjaT)U=i4wfh8D(9=}O36j5TOO%0l&^HN-u;rNv}kc;Z*)cvD`=@|X+_;zi<64)toBQ9L4S%F)62&U>uo*6(bmf}NTShe3GF6xE<GcSlXPLZ884PLsaG9(X;i(nbdw95rB<7Cvo@s<<S7_IPH~DgZE<)>1v7t_zS!GSF$dTMUaVXM%ac1&PHiW#^d|M%^pLV4x1>X)jnDj0{>DYYVeg~_w2O)?YP`w3%k|Y=?z$UYs^-PgrY{P-F&#7M&|Bcz>4WsR5~ECk6a{oWmtL!Meb^G|+DUm19qc154?zm7+mF#F=om!%J0a^XWnCj*Cv3I!g?3QpVCb==FLZdL-j#0FRUk}o2&g5U&YP5h42L#MIzj^V7b~02hF+IW0k>ly&mp%!{`G;aik3%_i5^EgC_Tj))Q#tC>axmI>5|I-FGyHqA7q1Nk@UPv!p^2)FLZ46Ov6_SRie|=tpFt}*B2F(q`34#t87+di*&Ru5=x!6h$fLa>E)))A`a<T&n{#0d||`pWz2k+8Igc;fb=e<5So-@vVfR5Elb+1CKG^wF%#=)L!~^l`Aa~(1AVqaL7my<vTEPti(x;cc%?GEgVC!FJA67XGFU+#i(XFPVAY;5NYW(I&=-{<#iA!F946G|O_u1c?}@mc^32++r3#t!=(=8VB4iZwp!XKEWtkgO=N?MO^H4eKVfBf(cS$nHe@I;xBkFo?5%sw`-96pqL@W_~gpg4Fg}&X9f)ox>Dt-+oR$+Fvm6e9BJf}<uTO*~EYP^cb4mG?-DRZfjBK?cL8*&NcWxrb5WI2nJkM4%{ZXlwbRX!dsS<28N$|;rgmS2MVqFvPS(23N+?)V`2Jf`xGrGyfef!NR@>bpCtCGmQ6w=JU1NRwmlQ<}BSNmfN49oa+$OcC9Q9r)Qj%y#e#CRFF=)?}E29CMjUow;0IEtM-Q*H_v=7cI*42kViH`&!r^ONhf$@cLy9Mn}Sg6_<sVnb9F$5;ioE&YMh!ZugQ<3NzRO>+H(C?%p^(mUL!2<}TtLOG8NK=o9py_<K^1ijkz;(nUroyCjXF%UwBeI%E<`JHt<c&-l}&5=F<CGpuitH43E8WO(#-QuW2ou8z76g7#51SHHQv!%>qx&unwmE7>wi!)hO;ucS!El2HdbqipmTa<wJ1Wt8JDw<3B^^?#)WWPlax=udUbr93Q-yq=#KJQ+J#eRQ>so%XRsK34W3;~T_ROB*%>d#VGh|9g|Jz%7&tl1YztNThjWRVueh``8gtUXzZOtg*hVwJoCDD85lXXWABn@Ufw80(w)Mhl)9H24ua;H0rp-@*Wgv$To5t71}tiPdwc9tX9@3Tp8GK?PEjnDb}$k+2X6EjU6FxrlY9qaw9g7n?p*o5wBUug2&j%94MAAw|7-xmA*>PD&5p|?%(ON0`8)8ubdL8m$-%0&f`!ZtBg{%h*;f#A7K%z9$=6rLv8b%x;wfFQi?jgS39p5C6h_erB*nq6G%8EL$8h0>3y}b;bKTDNaaQ1_Aa#~rK{JZP*A?MB}~|rLV84oF|v&kSb9YJuOF4NC6qSM){^}zr-T>{o_i1yNWof0J*BjVoR?)JWfG%Kf%bz34JaVJW|=LS-HV8NMBOLFo7#WfWNDJS94cKq`4cjJua@fU$;QcN7}{j{&Yf5UC3HC6aVVgtK5(4b{BG^srM8E})|pRrOQ_G0(;1vw*mtqvz#*aTn#3R75|Rs)T1>u~`gD1v@|L9DBbV{$(>#%y{h3l#Fqw7Mw3w@95#`1Nl1=)_c(s%uQ-7hY<w$Ft;~#CS&Y={$<F|+bmv5{Sr30+9oN$t`B2uhYsTNV%k@s6A&(0m-QFI<#=SMbEYMG}cWUQz^n~iACrgO;kJv^uMA$O_BpC9ciqLjE)ssxojE<UL;BBWsKPQF^IXfBeKImqAWg*V+JJPpNn2d~T`c6;zenxDAg7VccvhNpp#8n?bV9BHS)TiV+=WX3o`SPZUdw2axB34hn#4(SQ$M1tCIv)R>~kwdEEcX8m-A+J4Uyia*2Nh0di$%1#H5Lu5GD;oxbqfH2o&u}N436=)YEg@=B=Sf;V6K33XiR6wt1Ec_Q2`HzkOCJ>-o-|evn}`FDFBJKF`drJW=s<4en68dfc8znUy)Al5l1OZ~d>uIi7XgEmE-6&0H@sfG4Ml52`Z`a#Qm>Xu(UT-8RZmESTJ<y}I!O3LNnr`~Kzh?MYW8ULQAF95_)6NW?md*lrwx;!<jGdtgmQBLPz(N7w~wkIO#f2c-fT1Ugz^Ki1LSIWP&xJoR!VyQC>I&!#uFV#y|+WhQoF-f=_LnZp7Ls`)R^v%l%tL=dtGJeiO>tB9c{a1OvQ0%#z)JXB+e*qY-Z_pZ)(+Xm&%?df62(^p>AY%rE8=!En{~4F`A8roGtBX#w<#C@#zg)(<|QJHFhUBAi+i)5wMF2L-j*EWwIdc)l!@{X;eM44wZCCMJut1(m#4!+eMihU23WvC{;9<eI}9$IwxB~c~g43N5f^>2Q}G;!F?Q9p|y)rQ_YE4d#s%op}WD;`Tfbwo;i4_!Ha;7b8*h}p%yTrNXW-3P&(p^l_bK*7U9T{y(-#y9X1S=jzF_>6YY%XE$Atd4S*a7{h$%z>OyS^wX2Gl3R-QuZ&D4M{^ay5zBwuKgE*?=70s|?nmKPz=Q;H%qE&8_Q+bIr6fOjc3sS(cz|xCF9tNfoFNC!5^OX$+CBLm=NCJuBLYfjIPG!u7UngJxQxtLFw-S~Sz{;uAohF~MMT~@Hy~_*NRbC?MfN&GZXKTxkWC@*TmzHTW>fnh2hDtt(o|E#2UQ)2ra&c{v)P;<Njvw|}7DIXswWGJEv*ZK*bZH|JmpMC*g7G9_Oeq6<-bheeM2Ce*`lg|OGMPzS3Y|$g0_QOegzU05cM>v7LdA*)A|h&$(zAoAmq?*P*U=)%ms3cFzuPS#Pe(#BqG2x*>M`>naS~^sIw00zQW&Hs!$*9vlqgBzm(+{4th;hWP9u|0)I-coRD>Z{LfcGc+8K(%$5b*5x}#e}vfc`vs}4`ux<(FhntnA7hxDOEBmu9ZX$zvut|&|2FDulGwe`SRJ`{#a9bfFe%$AP5%c-OrZrZ84R_Kq}r`%sHZF({p;f(dK&gEmErIILpNs(Gcy(B_%Sw`I~mzU&D6M-EhPV<tH6j$bxbh5*W`E_!m?P4U;slBsrlX1HBkE6}tQRfFwfa}eT;{NKQEB>TNaFJ00;Eu|<)M=98my1yeASF?nit6c7eVcslL<U{@A~}mhdaIDCA5Mny&g53;?&vI*z+tAGB>kXEr$62zb_(}8hVq4WC`CeR#HoPw?)1$TQE{dYKy%a)T>cEi2RbhtRfotrB8k_2(Wt(MuIu3!rkQ@z+tJ;`qV-ox@uNDa=7XFMH%%9c^aKiZA5~fUjO8XvD=Sj&6exA=DNRf-+%h&Si!K6MU`73oIuTLtS$81#q!R5YxNUK!0I;@tLi^%a{zUbyCK8>kU*Zfzk>SAQ!nxJamyyt)DKv9B`Nc}^8F+`>?Yp8&3W>>sP_f!S)u-qbP)==_CJ8tyep<nn>Fc(D6az@gNQ0aLqWil2e+oun9#cJ&qBMVf092rlNn}U@M~!e^8d9#{(bL81XY`#mL2ehjEdfbAbSYudWU_8AmP+q(pE#L0IT#)CGW3!%IyBp}(Y{k=!rPJ0iigLE0?wjaEH`RPsO{A@HKL5VZld0Jar`IGkLb&m-d9TINV?E6Y6BaU)Pu4t52}Ng6gBBW>zOVx%Hyjfg~<5I&g<6TnUS>BZ$^E!R3{IIfq=Kx!%15Ce0mX5Qk1Nwj8wg-gUK(Kh9C;nGs)k?7YwJ2y2+JGru%WMQeP&Nw*Rl=6EKLDkvxkLQceIUBqi<K6=TZ_yP>V*pip%UIARwe<t@r%)S+VUvTGGO1Brx*zUkFU{Zc<s<io3(*;a`m#lYd6rGNQEx+!wpsR!zER!R|xq^bgJMZjCcM6pi~Da*PQYDO#b-SX&lgof*V2Zvnuk3u{mN_pr_OEuO=JyQ=VXHSZsN=UnE8#05WY57l_@B`JZogqI$+3NX9Nr=BM1w)n%-*X=ANTeaPWfAoUl&k|Lv$KoT*CZ@eY3ZmfVj!a~JFEp^-$i+?oFI~z(o1bvLU}-$-&sJtnM{m6I>jj-QwLs`MNzU2#>KhO2|Wlbwf}fVE}vl<<Xj{YX|I-&VU#dVf0aTeG7&kt4k~@RL<DCd6?0VbkqRgysYH?`4ZP46k)nar^>hK3h@^>8=5i4SDb~4*QGsLCp_m6%e$7F>UT&{GPJc}5D3gF4c+eLii6|o(Fzp?2-w1cyMg`^BY3>(G8@?@p;DOH)oZunxtz9KZv|-e%p-c^wK`&c^w>!WdPPDF=HfKqhfgi{&wL4+Rl2Rr#u#>t72Xzu>rxKSz0YRH*)ptl<d85F`StWxn2Vdry`tn7}fq{{V0g#ob?PiR$JguPMBg9uro6{eS6wk;rb)ti%*Q@Oh#g<WFj}|51Y9zAA?~r;&Es|fWr;O5OxQNoa%lMO=w|8yNcT8-KjAg>9T9;FPr@TR(5@R%+R7erO6?{v+XjPIE<o0z~4!phEbfmB)*kd{Oua*Ke$QZO)!{k^gHcf@{GO3CwB4sboDcWi+m5d6GbRmszlM>!C;vy3Y-j5YdVf30qB&qH=Q*MV<Mo1IxrPoJYMW=`~ejqVK(NiBp(NWs{A|hd)Obn2T0|_TgfZTefI`m?t&$0Gdmn}PK(GZU0dL}7h0jYuNmw=!sAl$l6hk|OdpNIBq%SX~nyDhhBA5%F47LO1%%<8rr{N}n)Sou^UA;}nkbKP8eb(TZ6U|iu8QWza2r&B~i?Iv4&=T;Mkky?7Sl$~=z1jl%w7^K_umQaw=6=e!2KODe^E1nerb&BO69STpit*3;5%y8)LnARm?0LI`*u9D?-!cDJ_nw|b-kF`?~&>QUF=}1=I$x_R-NH39)x=2anB*R4BuVTKazIE~(odU$mr42t*k6KZ*c7r8!wIY|EP4;jOjynl9>5Pge(eH8?08bR3n=K@-ONWhOZC*-MW`jC~J{7BSEU_eR=wv)Cqb;@GJ5lBEgGS1T?J$z#l7o{y1*VWYHRzTl;~__5*gcc0^>Qh7N#bO1FC+3a8C6Vmyn=I?GU|EolXa&9pof0vuPKjmT*nR?Gv(!oanv&}8536r8FitF=cbHF_h+~e;hU~^sDt7orVASG4gx@^k>tN_BzZrbvqsuM4=dB4D8^|w-dTIKwBaoToB#kLcJ@QmBQ#<Iy|0CvP1F?iygDC>x`*TjQz%bwkbB-j%2StxlQIX#y`*&E)45414$U}HK06V+q$wppJm5qi)2um2F*J|`I9gXJ6}Reg<5mQ)m)ZxgZ4?Ldd_B(=HoX0V%7%;Nc(j3hj$)IFA)hivFp`NxobCh_tZ1~Be%7hOZ5fGpv>k*Odn3`=NT|J&Mw+HES_7fVojt}@7fRTCYvDm9ZUmu09rG&TkeYD4E~9kPEErU0$w8YvUnmKIe(I!zQJrxORG7ij=)Fvx!#URDP&?9Ru!%ah5x+rRLpnXx3}jdLk-TYMIjTaTcm}Mx^aacUA|9yN%BgpeBk8|%DP_@s&jlG7vb_h_PN(D{qV7ikHW8?ZBV@$N$rE&7)5EK!c%f8eWX7mm+ey%QpgK}1WlSmcOT9w$yP|XnNqSe_6EE#J<Hk&)f($D_k9|@reh)&<U~#6YsMF%1&Z=tTwQCD=<lRP|6xUx)0sh9tI0C^1VaS!&!MZ8(yc2qZScz%M_+(`xB#>Rh*Oqm$4UOb%NL9!{Oc5!s)kJwH>k>&u`NLAjNKU$b$tcHybX!g2skV%|ozZbQh}klM8+uQPlK~168A}-)S3p6ES$R^t3Va{XF;^j_ctLGZKtLbXqbvH6G}Xwn@oK3q7;lhlZ^X0))Nvw2ENPa%9qK6%GNxKiYD6hGvQZ4@U=n9FTS{COIcT^pD|Xbssxxdr$UufpDM3vP!xvqzlF7`}QcF<4M*$2llwuPd-k;$r^-=QFKnA*WEhz)Ai_r)~wHL3J%I_p|)PX91@(B3Ok|+T4TmY5BXN~lkJO<K90SkZ*QQV;@CT%H!B$Gg=M6&!<?UHkW5k-w2vyf!r$=p~*JShCAD#a5lqckqDDJX?;iI2fc0gA!LRE%)b)xS%d%Qe+G(DS@n+Ppz_B@a6bz#_@_9{?eyga^SM-abB$_fcADP~?>#R>GDLEH&pr7q+f_1oV0!B8UndAfrjay*H3ddDd3tg#y0-#hZLuwaZco3fgdGuWE<r6!<t1Tn;xmGj}{`0;Hk){Ay{#zQ}{n%i#2IsC7jtZtDM^2<F2uUYT9PS_j7V0OKLeEo++|hgHfl;rHLEv@9XzawOD_tQx4t=>`Oc3E;OZr1WzID%Q!}^6_NrG8r@+luO(}-bRlOCNTnED*Yg@N--wr?vz(c^+r1ogn0Pa@h)sA<S!j{mW)AYD<&({B&8C)lNInSfhakdEo7m4qL6ev4G0lz&5qE~7HvLmut}7S^J%3<K}`p<ASZ?=tu{`9TywrSbkqQJEV;W?r|Ks6<c_)wRGefTUM-~vt{_bJCW4DZNT<m<mM70j4#kWn=^C;n5J~tagEuH2yoB?-%0LlnlhE0PWBI-^HQ1$w`EN7b!LH$m&4f-R*iMc}+$y3+1x<8-WYz;;w5y&1Wk&@S6?D^v+(FTk1;KZEwbXqo`k#)Ov{Ck-b<8L}5-l0;Y}N?2zy(v2;#FwRq=E1lq-ajP!)w+{?Ao|}DKqL)GLm;W$Yxl^2%0lHe~jNxNig{e6dzZc_JHD7G<1<szBKuLhmunAGChC?h`%vMFP8!fH7eH$8{VSMFW`#<*jQaMl2Zg}1iYjnV<2I5koW;9vSkFMEaiQGW?V7?8bk3o%q?*#o=E5d&dfprKcUP{{?dUmNS_GGjwt7%cqkwmy~t)tO_w5Hk+cT4{BRoK=~CPbx|@I(7n>+TFameScQgKwD??#R7lsA_{s>U2(eqh)^4v&HDQF*J^nF~xx;sv<ZaA2^D%n<&2~<vy5(}+kNJ6R^MuATtuZz-@{^}wn-cFu4^`xLwujfu{5V(IDjb1HnRIw{20VEj5m4IDJg679Okz)ra2D||@aZJ2ykV3ROMBRr}wv_nBz?Wo^^D>iIwi)Q_KJoo?FDXUKLF-e>L~*W>G0XEMCVE_eLv>077byph!$D6-BBh{bx`Rz@WYkDOy;v$Yu-HyH7r>vzkPS_?$33N_$3ualToyfz_8c)Yc`)f89gc^F<qXPmq$w`0OsfMGyHbNzfYPV#2?Ud?wG2~(Rygn{r8szq0XCHgR(9(lB|u8b5Hfwv9n)H-B@|D|tfUt!X^^0$j&cQMds3IPF=0@8nYZ#FshVll%CuitNLurx#2dX3TSihZa#sgA%Dt7;lAMu<HTPNmcZHQ~!c5S~1x2-j4j>Ef(<J?d1v38mA|u5|#T<=mA<d(=NhwFxF>+#Gt!%`wI^ltGs2^V2j}ERCZ9PRSpgquEjY8M(izw5(rJ9wOjB-RsGM7Mu#tJ6gq}H910QppB9w1OcR4kWf=^VTaRTG5<$!8mEdT@@I-&};0$3<X>veZ4h$#uI^?m>CAl3y;B=QvzUmebyxKq5cDcB7(ZmlG(W;Il$gk?$n9ESwgcgzfW0E*9FFR?W4M+;4ws<4e)Ee%wLH7V~QbOC2VXXCV%CK#57ufE&Db7b)GuMZNkUl97~TI9+DdkFEKiE^Rms(lq*8QvG(mDcI>qV87R;G_Is1!~cB4fRvLX$1BL^UMfbB+RK3Lb{qmu-|v`Uo-|&bVqz`H!OzWXbUV<gsmo=nj_pe{yJEmE#*##u6a*<3DY@YeR{T_|TcO!Wf4Wo=fljB~B}!jznjLf*8OI>sXI2j}?6%<nj6$;uC;`cURxog<;{)1~+_3bXVw5Fl(y`^$0I3AbI|X&JRo}BtEiZvu6o{n5g32bzg6Blm`y%L|t3OIPYBwub6KR2pg#;m&tUaGD1$2u4hO3{QqZrarG)~YiF+s!bH;`h%$uZDLX9J%mkRFlJ-cr)>Lr3u)&wELkX}*#H1-V`7l2JG0a7(0R<X*{dq;!;$ACw-QMv6jPck%9qeHwDWgv&uPraSdi+BvlY;nmUxu2E9z{L#;bj3m`!a+EjLo2STcXtW&0;Xpj%*U~EDX!Mja=<9W;fO%{MR1D-&lmn`Qu*#p<NJ|e&Q<hU#Uat^9TAvq8)&W$A0yUt+7dbPpg3Jx2fN_)Fz87DzqdP`?v6KdM`AFGdUBhWm>z#r!m5n8*b^%0;%lyH%0f<+3HMym2FG-F&hiAy-60_Es-lm~~{@eMy)RqQ@XR06ob}>1+%-&8K*L8spYQ5L~hXFwP0A@*ZaFUp~-Eujc1ds96QYnba{%7)V21J8$7!w#lH^_V{@FrX>`b4JIM%#e#8pO<7PJJsCu_li$@9p&IaV7?^hD`JPL|P;y--Le{Us+LbRO>A#s$3t4^}9BmQ?tW8)7dyMcd}j}Z6}zS29y^|8%eL$d}^aRgcQcmQmK&J-sZeawUDG~4vL-ricYoZk`r~NQ%(|Ry%d1WTTXzNGmR67vfF%M!Vt=NGg-hJ)o*|a6QVgx*WpUe{nZQ34ETh!vNnbQCj)3rl>I47mGNp_i+DvKvvQv4L-|2Ew5b^mIrTjr8OcLp;+&=dCaGxHQaHC#&LVlvrs#d;?sy#MH$e2XV?xUhlC|Wt{7xx5@}6kdA)DZ23!k3Nm{7_^Kr(M2#o$HJwXBNG-7O|hMZV{&m5p8il5qN_<Q}0mQ!|)E8J^{J>dqIyEbOR7Y@Gx(qzR{l?upQi){Y&k>|;#dOvIoJ;tVJ8S0tGn+X>|vlSu#%NP3Yc1*mSz1uEkrqk;(OvPP0v-lwdqEFb?Z&!t)`K3z%ht_#$0qC0w~;b3yvfOy6gCE?pu-w2&0WRz{CrNf~6Q%1LE(I*s78ELn5+)vj%>E)@(Tn{DnY8hpMskvp!vgtOY<Z<Qnaynf~deE{&Bz^`1UW|29q(!0=XzDQZ>BUNV*dkYo>RLc0RH{R|pUt?GnVmT}nm|`%kZnWA(z26ubnWoj-Z6P@4*PO^I55$f2a<zr?wcO5J4~9OI(#Bxv#9~~px{&^<w%)bn}IDNp%7_D6px&N>}n`SfP!&;wURIckgHpgef=xFQDzH?k1QLpikzct@F+e%x3^N$R%Ek9)C&MWV*9azXG%a0ae(+y5nS5{q&d*sV<9WOM(znxv=-?6B8LrNlCoiW2uX#IvK*nY%#!??J0^`@vrd!mi<OP0e*n-xbKQ;(Ix-N<qN?%?9%@3*9;AA@dp+C?Z#Xp2$WzK_h;?;*wWNf>1R?d;@~qnu5|r`4F|v@(epAusk0_ief~7t}KgjIa*g_`isUTA2WJu?v$+!?<FyXGhSc$pOi-+pOS>}UgVmkgQ&F2*b(qp0P7dQdl7ws!GUK#^*te1@Ry5PNaEXO{qr120^j@SUNG@q6RR?urRLqH4-<TH~2J-Csa4?IY|xnGW}c7r&EX1D|aMFjFszIuSOc)3!MPxpo70et}r5kmx;Qx>@eDO{>aq=Cb*Qvpbvz_TS;zi>)PX-59So&3U*vVyHd%6ga>B~2-%Iw~2Pc_~41L<+|4leqP}lYzD>w4Gk0RCt=*%>imgQfgDkP1COm#{b1i;5_7=Ar1(>gY@bi233O{oO=n?%4|6<GvbQS>*!@gaJk!gl~07NKPvdd<pt<G%-SPKJzW23K!~FQJ$Vc43Vuo|$~+WEsiPxzgNUfj*3>0s1Zzc706rca*ii2FkrEzv&aYNt{0`uM_}lgafS>{BFTA%2WrXYqyGD)6(GXaAu;u{c&u9zKe+P7UCq=e|>=Zc>5@;xp`6;AsKDv22?3r>Aot3l?DZc3V4cf0p_g1}FkeFQQkJYa9P3QFPM+LQ$U#%lTrh%F|Na`FV@GGP)8e}8$w7RKIZ^|hnXXV}R3clH`Kc|S%9Tx*<Z|fZ}LVhRo&g7X80N_OLMwoCjnll4_y`G((f=(ALnuU~t-t<L8dT?5&8LP-O2N_xDAd&kZ<N9J90bEVsI}Z0m9CvCJWFJwa(^n^{&z|+!lx@}OxEf6kSIQO<h5~GaRO^Pg08+5`LR&Eo^K@E>1~dPjyT~g|?&R%oY!Ikb(lr#~Wh5>#x?vL3Op?vLY4h(cxE-zn$5$&`ux5>TRmq6V>=czQB^coOrlXPa=mwJX#V3ziG61_Yj`OV6U^<rLj?G$9LfllRBTutJnLBCP2FX1=jRXk`poe5chMtw)bRED>8+*J+nG9W_H0tbCa+#7|coOz0zFJA{%Je8S2iqH=GJ+2aSg<T~5nac@Ez7g?PW^>S0F*YT0a1ic9TtpPUyeENbY|g*D#kH&4wMq6L|7sM6E9$+!4KTpoj!3an!<TG&=(O23lzRLU<dV<sB9CnaIC0oyja;l>1Mi52<SOHju^qp=$~yCgn<d@ozu-k(PF9)r%MxUaNDH|^aRu0d<UAogyc4?>=$ULJnaXdrsm#0OW93JmFO!jbXd_GOy?-XjxseDDT^5%aBvVR0RYWzz0pU<@M<MMLdaEg9O<G&%0Pc|cE{i`;S7=fPLi3FE`9iDcLWn~nE&FG(ycwvz?~S-0Yu;Fr35p<w2DLvgPybLxl*ek^-W3<n1P87@F5E16YVM_kw*rTQ!pR8$8?boQgb9dU#)aXl)xPbqf;N*5h<60fOMf6sRaR<)j4#Pk9MI9pkU}ypPEOmfnH;UGjBt(20+;E1k-`Xxm?%eqH~)JZZ1w9@O_c_8G2Nbn+GL7Sa~k?1F||W4rDjvGU!O&ZE?Wn$R591Ns9wGg+!maK~jmH8bd1kF~JOuYjn~Nj;o85w3@;V0#LG_dM71!lXo~lwsy5IIICgisanZDIfNX%*$bNB=ps#z17E;H^V~$6REV!(7WHm28f5cj3sO~F<zkVZwbNk#d?i^?r^*_S$a~dF812YLVSX&73<n1hP(c~fd7-_vK!|hJd;)b}tO35n@>~ZAGh;@EficH?TF4vSKkX};>=nS~P-&ZTuXI5cy6tMX9M?=A`~iW79;Zl1I;WZ$?0CPgR#GVkv(BaO&blTN%Ndmcz+?IXbop{prT?O<D=^rByaBIwx^(azft{jGzAcEcna(R^_c^Scq^q~*l@my?KPLWH31e=GVW<R>kb|n`;tOQ-d4gw<*g!(+alu!hSY&fwEHtdXxr2lrQDo%YQ%1u*n7}EL8L$=dE~0Y+iRPNn+<1S2K<Ff*TgfMWAJ;JvPWXOp-g$^U+{&A0YLd#K!EMCxid8G1ze?9W-CN~iBzY=mWT@0gO{YSN4rfd%9-6TiE1Nr%IY*H~E;7niDoDg8($H%ehw*`QD-#YOBhV*&7#h|dGLmrzdfs=2juKXmGqD81@tFOCI^)!MaObC#5vLB68bC)nhWtnrcIn=~NSO`lm*1cV(v`&RF#SuXUesT$)K+6s4j%zWAe8a|^#QP5`B3O86Le#eC>c9TuuSZl0dt(kl|qXyO(to(gbXlWI<rJzT=!P;Kb+9j98w6|Y4*md=utBe0*q9IyzU}og;~b5mGm3vn7ir{@l65py;#}M+yv|xIj0>V+2OJ@o$5Ft<4_ktr&!5RPOnZk8(;xImWP=HkOouKeG>-}BOX2u0NiPfT33=%_E_);tPu?y=n|L=X>W6k;g;@eey2c~!Y9$fFSV2D8YGjGr~<!SiPr%_88~Ej_A>&pz#Sy1bIdo*(*vNGFndHaspMfsfGr|y03e^z9R`~vqd7m25j+|p<MhwtVan;1P}Y+bc)aA_gA%*ZZmnI-@D)pSK7>TuV_zGEcd4CZQdsc$1UuoyLZbr)mD^0I2pQo*Hb5U3Tba***#<8e!pyo+Cc{CZHrvK2Bw1Mbes>wEeI@fu>5%5Y8J|8~fKkHFF^>sBYIMaDJQTbK!!QPz49X?FQ{4tZr>wizbd~R80%`HK+pCqcJ+?r}6iu&4h)C-cVAkFD7hRZ&w7>+{CV|_T$rlR3Z6uFwwgMSR?QmQ&l2GH->FG0VAj_!(oeNfPtP5|@ZNhm<Za0N(K_`y(uOnWBG>DqxP|dv%+nuCn1r7>(;fs}ZY9{(N>!dG=CPoHAgzGm?6KENbIuqMwDI4dMv`RT+5YP=SY2xZ`DFGI8LTM-_m%9pN^j^}d(RMRaHlqUB5lDNRPj02{L>~V`$`+bI%p@K4Ga3B5TX88>BomfCT?uVq9wrXp#2_Vx8F~dBqcP1S&3fe1I@!iN_tG{O;-h$abpXmqayio4=T)0zPma?{ckqg}FxJgfFVha=7oipLz)fhKAr={A`eaQnVnS(4WIKTF-}HUp+$eD+8hBr=tN<P(M+yI%w@y}ocJ*%T!0K|Ok;?^?FoAyb;s#Q&P*I+`r_oKPn+O1MiwU4)mOKg7`m#!J<iL_p`i#ouKu^;~vOSxZ0)!RTfNqzIm~_#@E+?gmiw+~G)GdX&LBZk0N+~{KXg4xMcNGJ+;wk9)Ny~?(bgKJcZUZA<VH(PyI2$AeTg*yUQBQ$SyjAR&93`EOSECxXA4*GfoEkjBvXkwofC7eS8fq8lY4z$a7n3Ty3#TA03f%g*IZ)D7fP~N2LHM1mk~#_Vn~9z|S_hr9b-rkFwKUF2f8;cyW}-jR%is)7<s;JqCbZmXfS4i%06nv+YjU#fy#nZ__cX0AhhDb^SU>tv%_~I0sDr@1$Y}Zt7ziOZO=JZJz`}$_9$zi2V<MzNb0Nb*Lhx56m~h!^dSr&ZJL=(6!w&j0%$j!|*NO{BFI?Xod$<U*k*7y<a_V0@TQjT&XL_a5iU^_u)=8u^!J@o92q=Cgr5+dr70$MFaXKK3AtQxeJmRaByhY>2s*zHNNyixEEfhbNl-)TUN=0QV*z8&x>H0J2N-nyBin`+kE+OR~I8g4coT5c<D}M^obtiQTrv}C<o&+ID79QrpDc5?BwThWn2k9r-AGhaHIw#%tVj*4ba?;5m^yz_srHe52&Dx9OJ&5N?Hd2{hZk=NMLXSo(r-&eEQURNE?0p1w%V$t)BjsI}rG$|WwY548O(zsh>tarJI)G;zFx(Fp8I76JwRCY?GC~dn^PS;SuNLxh<%=_{g9;=vLbjaEW*HXVIN^~{6v8m}c_kx>8LBhkbkJ${BPJ*$e3zx9*UdtbG)^b*gY>05-!x%_cBm$QDLt*wq!5jF`kKjDU4&F*l?;4@zd}ON%S3t`7KHR_A?PKypa;W)EhMF*0QWEhG^UUp*t~%auH}P=>=+(!kmq5|VD$8A?lNYljANEOA@ri!QP+*cO#ovtcM{+!6@zT%YxrWomf{Gzet!{i7-SbFm#~B+YH}gHf~tMBvI1t35rhVKIJpc1DPqH}yE;~qo#{TQE&GZnVd68r@regL%nzc=ch|o4C8m4BZ-%@BK<WKNZ&b=}OT8wOoKfTF)<dm)m_!pVtHg(x#1pi6&}%J!@TPkOB#Q(YjN|pyO3;U3o+qAna)R6LA4b$HH@apA!0<RVM<VpCV_=EbA5xNFQLug2ZB9!{^tWLMI5}WS=i`#cgMDea?Ih_TTXD85aWAQ4V@k<L9PJ{e$u^;9s~{1~><dK)_b!`F`{_!Wh%&e}YeFj@ondI3qZu96ZdQ_%<Z$3WbPfa#1x~T9;;DWTgwOguIEA<RY2@h~5P&Zs!dY@2?hRYW?FigxnBxh<X%b@80{oM3d|c%096(T-&^X*8=~q{H;a14$)k>KU{LUKP9cfx446|!6^ixjn4r`-b2H_5^h2_*G2d4nca?5Gj2B_<yk6Uu~VM!qFb9HIqFCpe=cK4^4lp4TygC9a$nv?*=PXpm2M-bsgCa$1QS|=W=m^hpAJMhg{E6FV6FiC34<o(I+&5>@EtFw%WULUL&hD7(doe-qKNKY;P6jKL?s`8x{t0g8Oy<r!qSvc?0l9Qx_E@E3nLdeN+v+<vJ+6Yd?fxK?_<=6&n0O+ae=svcH$*-m-O=gFEx)SP5d4)tA`Ix-T^cj;OW3bW^GgBGS`-WL!b~z-;Oc}01u(iab`IV{~y>Uy-NLPoqlRoVcRxX52QEIQCzb$1XWeHlLcuTZ7#6NssDc8+JmBWO<L;(poNpNwW=>&PX(yT5sOTFOJ!Q%mprcSf}F*#_V@q~MAMGLisJIBYxQwLMdj+*6i2?_Zdg27i*7}~V<vL$zcmg!xgpO`oR0%g+m2KRN+ttuV-kkTMB_|Pfs^c!4i(jl;u=q<inSPh~xqK^8yg>)J?;hxQS^MK@;ia{X(Jrm8OtSg9pYAC0a(OEM^?_>Kvu>Lw<M;c5X$LhhS?;v{u$Rr*h!3-XZs#-Cx9i1!QTTuCz&725M&#<7gzQSY<Cl85s)0aG5SgCn)Xe(8;Nf}EG35p&vkgL*WwsiREj+q{3<PNb}rtkT*Cy=s8eA$x`meq{YE>2{C+`x5UE44OfJNYAXvv9!Z|9~@xeiJZ8w1AX!k<sO8RNO(=k6s+mU5epxvHd<>X^1Gd0_SCt8@oD8Z~<^&%9y~!MoK`-lcn-FU_*l3<&=>j4biiL76+J-u+A@&yXtPlSyIXh=x~R5Eddkr!Yb&*Y!0i;<Q8lPgc=Vqq1hvNBzlm^ES+YDOl_z7^*mk2I3>6I6A{nZ4a~wQKwC%$&5)$Sq`>PrR`%|g?WfXn$73HO&L15=y4XubXmI#Yd7kbSTz3+39|k;431N1xWIg-hnFNpA=b)>Ltt71I`11f{<*S3OHC;5FqDOna&V3WfaV)jm+hoIn1%PmOp>04_0ps_`iU~!BG>UAQ1B)r)F^pEx-rsv0TkTZJ3`mGXl+u=|1a@oC{RAp@eS0VBH|8Ba)NZQ_FJ&nGB4u+jg$y16v`!Y~o&pfnx8lXh>hTK9;DO<G4g$bOfg8w8kW<bM)qr-FNic?jaM*}Rbdg6p<%9-auL(DIiHSFF#21<vbpJ`t0#^t3Vv30sCnm;&^urzuT)?R%(0zzGN%uj&k|evy-&RZl#%rQUUahQdL!sq#6=!HGhagPGA55O59OgpEc%jJ|^_bm-N`&1vxR-hsZ<<j;ZP$Nr4ZN8J_bso7JxJ@+zQS3^b?`(<M$=vd9~tCVm7xw`OD;kZ_u?3mv_i;6Ekp*vXx;H@CArb+LIv$9YZ|}^PQAO2VoE6`Zx=out5J^c!!b>8GHxji371I%PO_yOZ6Zvd7)LcvV78~Vklx9wPru<lQr=Sx?x&UK4>s&Dm2$nvDK!k=A1VU2ILJUFsQCh%y}w#%8fYyd_Fj|+<nQHbIU};o$N&zAwZrl7p@m~I@*?BLgB<vzocsZC@3JBkO`UKqbT;y_%Z@;O6s&ZoobcR$6*Q%4qJTj+U832)T(+~&l?XL2m?F6Z8D1R5d(!ZnUamxA72300BJ#a4Y`c?s^<>T|W&yYzx~Vf}j}`c=4yk6^t`AehXbcYFF`+Af{cu{@trv7puBaq4X)k+T41Y(Lk^Y-YN;N9oU1&ST>|4Ibi7qzkfJ}(wqzueofj={?`m==|&R3m#W6YSv)cGNYd#v7di&+V+A@dkvb>x9jr}GQKic?Icg~$<sMs|xy(ZxE5>(D#^X^NT1sno>{LZ3E_9EMp=s%j=V(b7Tqcrlozn)%1BW|FE^NK0ocn*^KktCdayK#&H(&%7#bXr=O9HkB?&7+Nz+^;EIE8cj2l@g!~I7F5c_7{z<vLkK9GT9fHTftyX!r!#BlAjvy6#{@aPNMAKk34nx&v?V#;FD8XE>vosU{FKP-J$M8olo1{#@pK_zeZnADYvvdPie4U91F&wx;wauQHYFB!olxJSBaGB+I&bG3ua*F#=6ySJoxu+6%#0^Hmj8om$-KiEtuvg~(*tU;$DU@ghw-eGC@x}BOzh@V#>Prc*d@^|@Nl%3E2-kbWeiVbnHHqGRS3$b$r5<rwjG{jIM`<cwmfI2XU~7TY$w`WmrTL*YKn1Vhdrz85fsOY0>Dv%`b{fb51yjJv1Mj3o&6&*X;c~%%D8+H6y^eGZOm;WK^c3bOTjdQYM%F}D`CBcv6lEF&GZ^nRbfuvWC+z1G{8g}Cnv;g)~ACPwhp}IQ&7qV(p2=4?sHq~90`)i9j3VG`hZea_j1gT3=$|G2nJRZ(^8^a=$XNLyCl0P^qbcL8^l;sX+3DUTp?-s)k^cr$TcD-SZ8Xhv820&<(PG=!}veOSqDn?#zPD@H!wGLXS<+<r%UJ{8CJF;0_<SS8@<?o;+G<Fg*ZYCU^xjKJvbSJrYT14M{)Y}E!PbtP&U6IwXGYRPqBOQ7(+GUtCf{(8*GmNf98Dzok)^1`jxGd3+1}Uqd>iEloB{1bdiknxuN)30GZIZdqW|3w47$&-`;yl9x}8FMZtJ8vPhlWs@W3heV;95m-1J1H=3O~@Z2|sGU9QOxq)ir)k>y}Mly-CX*!xXnj0)DoNqgYnhoRcz;>agMNB{<72%~Q&kM>p9Li+e^1YwDeQP3+gmCL@2JlGd<&HOD5e^<^LU43*YKL(=jqvzV93r0vYOHV_q+cT)a_F&UF8cOrWhFmFuABsJTFgPT6b1L>r+KPybsdK^69XY<Hru1A?oWkV9YCgw3-79i)ewEBG$!m*8(bqVZPBFA{xmYE&on0vBxk>9#66CSzJ<?a7;Ke6&&UCUa2X4Fm_y{nO7k}3=D>zC^QU!~9~WnVZn7z6C0SR<-N#9B16&XW09t>}>&cu@Iz-?_T0tlLZjOLp{WdrIKDS_DeDJ5Bh4G!;1Y{P}V>3V|;{aXulQ}dHYSZHbLl6ulR6(iwPR7}%D`7(dNXg8hFuxPm!9*x@$<80;W5%x$N8@n6Cmqc&M$pfDIz)sdfP(Ir)g>tC65r3@e(N%pzNJ}4ieM}~#0-IE$0PdQ>0Lp<Q?P>d>fn0*uG?1Of|Q_8^b(d)UhU!4N+XL{;qk?mvq0ZbASfuS9EO!-#0YjEIe*D4Bj$-DrV*U?vyiACbR?CXDe+xF%)|Mgyij>Oud^;w%IKmxq%3SJT!&mP&y=H4<RWG9*b-WJ7*!gh2@3<BDj<PZE2$db03A4oGjmbx#xn6${%~B{b^yQXprnQo%%s1PdjU#0y?BWI=<j(nmeMC!oo=B*!oaK5Z2>7HK>A~1Hb@G>^^Y?YYMq>p)Gmtab=TcQ3W{P-z)4KH3xTqlOXB6iI_C)+GQS|nQ~+>JV3}Q|BAPedC(h)_;T}$HE0qjp`Svvn?!Y`k*-58ZTt+1WwL&!bbZFP#ht{M4eNi>aZJ}YQ(Hu~HXz5Q&q&xv*UAA+0U`B(;7y!21(|b2NO=IkCd9~6!$JNv3z%=t&#dO$4aLs{j^poHSX4pX707^cxnI4#0|0$*0{)^nOxHAdaGxsbdSU+G^OH7w+AWfZOda9!Q4Y~+0fxaK6fhWCu5tAlRx;Su&bZUc<`5-hPN9nItnnEB%kAi3u^ia$^;j=*cYRYMjT#97Ur)F!U=e1A(PJ5~dx#|P&!c+{)xy_S_l<|Ps6PWYf&J68^2bml$wSZUQy3zm+C1B(G-~=_0^C2hm93brLhkd?K?j}J6sDT6dVr8X*Lih|gZSq_wV89KMC4|IriWwOcmJM<XZp(G|du3cpOpjB`+^fKB-PV_nFoEzm7u8|F77|`hXs)M_pyo-Yqq49lRJH?r!j!PDo7p@eBr&ZWHnza<X1FZqRq~X8rz^S0?hbG83i|*R7m?wN1N)tqltfUp^N=?|O6LuzGSf4DN(oQO!PJIvY$*q=4Hnd({M{%7LX|e#l_>%R;ChCf28hRCKa?3RevDcl9y~o*U{|Q8L)0tz?xcF6>--K?p01?zE>^ScP-L{A$-^nSZOlxKXAe|_{hGdQcvSV?Jy4j6&lXerO+uXcOl!O^yb#rnWbyz=ebb7tA$kBc*cJg7Q3d&@45D%JA2Lq_Mz>2pt>0{8H~}C{y^dfOT81>y(%xPytR8oXp@%bko>}Vq)Q>=Y?%ia>JDdl@V<3j&WdP0HTBMg#Oh`0zG(9L~i5b#5wnl<(b21w!!c4dUm+B+aS;dt%{UIf~;ZDZmx|LMo=)~2dp>HfCiZ3%%D!f`rR(g2qCWXDJfTUDGX9v{BYzU#xkO-VS3sCuS_oQR<l@UXoQc8mW{RU^UG!R-Vtab<TW#PP>!P&VwjaY!8o4&OS0a5ml0+cS7g=G3g&;nqW)@H^aNziaFCZvSM#2njCS5{{l7*Lerg6POjV8`)y;qsk=nqg}gU7xs|(8hs$1aO+eDQLDbIm}eJ1Wn1dk)AoonU3uSk@g6Ev2B3ya~~O5X7dB;J;=b!>P^3h2|Bi!sEh1<WOOe+J$Z%dL;G~4f*NYO41wUNS|+`7LPV-T+nzu(MSCE~1-27vB8V;s7sg{d8D^~~z;G@%sF({DF&SwL4}RHVnp!j+i&BjOeRegsy0}QJ?Y=9M3tZ$J&`J$FRu88*<f~B30#Wbg)~5?eWa@xcGyJXyliA^v1o=@~BYU87oSmM@y`9MGw}{a{#k5<)>=^xTiOFC_s~{Q#49uKj2AX_i*H&!|J8HI+&b+Fd2NpCB-1YT-h-m`eK-~rqoyCNal_R?6Ks;Mni{YhI1o<ip#qf4$fMEMw+X?S9Q_RT)B4z<BLzUPC^(`j=S#6o;jI8=t0m<r9YqWfeo60M=au*(qGgJy-XUDPRX_?Bvp5v%rwv&z;$|6Q~N$YspAow-LLYe~P>B{Q<6lMYkYN?hAUkRKAgYp&`2GzbwT~l^MVnQ1S>oh7lpYDNq*(JaoCwz%XwcJYR@Q{MBxC;)m6=TJxTf`6yj~?ud@j^}mfRdAc`t?Q5Lfs4dLmAa_LL&&df0`tYPgfE(d*Ei-my-rB`h$(AnnFU&zzoGik(ajXX5xQ<!S%eGlpsOz+Iv`1QY^NRBxL%rs0%<ddZOnXM?!9bC<5RGGmCjtkP~z!VJPE8$mWql^c^vR^>Ahya_mS$(W`~vaww-j+_)WBYz1+se4$X5chJd0%BjZ%WQuObNYl!QpKxAIW<n@lpcdD*;T+h=8|^yH%*4Ce>0|;>z$rLADdhWI2QZ1=pc@sxK}zX*Xblc)^yxHA<w21gLqXX*v{x%DqaPSW-7zWK=i&Sw4rc1|5c-_s03Nu}BFbpSZXyy;Y^R{mYU*~~VIP}<R<5C&10XU&b(JkBn4;RYpq{A++y{ETNf%h4zwL7#_(Rg@WG~4l&!b0D3^<CCGx=&|W#SV(BS3a%Sq|e5GbfX8>lS`+bs*gaCrnGLo{)K!466$dN#&_BY}6T`C8-COAt}A#n47JKXl=k#n-`W?*n=717MN2IpyNs`=@)dg2Tq$kDFVYhvU{aa3<l6YUq*caN>LNRrB(I=2Yt9k)>i5ZDpk@EdD7)J*vq)GPKQ6DPCOE3YTfmcx0IkmqooWOxouNq1KL@YJD@g)&w(8Nd{)3CAT@dZizmb=Em%-MXfQzKNU~tWUO4l%oljR*b_3pe3axb^574WWtV!BV#Y7DjJUo;6JhH2s&(O@D=R7O$40O78W<Bg+tggf%9Ih0^x5r2BVQh90kcMu-hu~=P){{vF^GQYT(mG85kbLDy(n!q6$PstY9kT@TbRj2<{<t&;B3MlFd(q>JkD9lb1u7wi-GIDf_*(beL|nW_G)97<Vzn96vc#lybavD+4>heWTjUJr$-b!fz()=zbg8ggP=<%uxH-fB5R^E&ctj(I8%g1pg5(amn^!9<i968ej53;WM51M?%yQc?QLmXKJ(G4QLwg=HA~K#Z;-H?93hVGf?s3yoxS_)_7kE`dzMaVoOag@xa7n5dw8LYQIVE#=15NvIb~13}A}Oz;0}cZeCrM2RVpXKwNUv5>U?IceOhR8c{!UjOu6{Ze^;z_A2wKGfiP7Q;d8x->qsUK5Az&V!II&^xmB~BSL(F*rFBi9@j0bLET~T<Y3CM`F&FWa@fFj?<YlKa^34>6LuUBw#Fh8E=Qa4bKS1Vx|JJ>m>qQdIJ=O|ILPi6b@9ymlT5K;>6DAUr(o*qtZpMsjOw^btW6|Jz5KA$AYnTF8w78HjalCvqPL9;QM1jeQiRsj<jlMK2T(KOqWJ1};}MIKbhE0ieOsxe-zte#U0H!w5xR*MZJo2aZt3_c|_&Tl}VPHv!@Am2p(f$p4BQVMf2It))<*#zJ_KD$spP-NRSG<VNrSjmNMG)g2BB44YRLD0O$sMsrU6`ky!Aw@EQ=}QORGe~T5y5;3cT_OfvQaq~*BS56cG~sF<uXT3K+XxpI_JTG^=(M0{q%?B6f4=QfzI!lR@IYg0BfQcEIBCnteNNlO1+BtB>FCK}rk@++_=$L=JucoK*$MM|LN%#B#5Ex*BK6l>d9{+EH0%>|n2rUGo;}B4LsttFazax^pDXUW`y>n30M=aic(<H~M+sFl-;&g0$i1V?IEkwn){-Xp4t7XNtn4?88#g)VV5cPil=7wPI8i2`qa>v{J4k{O_R+^jrsu^%Mv){VK&TDP{e+<c`n?3ydR<Xa1Aq=_O?9Irn9;}xcKU6%oMfp7g!gH&&P&%E6Jmw$948gr5|c{j5!MDnLu<0y4lvksf^?gcQUV)%5t9=2;d*=qML7)VW%3H;=JIMKeJ~JHbVB?I`SQzv<O%MT=@h}`)!k(T5fgocfkY7dH%`azFtZmv=Vn;}#?+mlCdltbk}WCiuE9$^r35zXj$sIX{IDRUApc<4%|*(Nj{#J!kRC!xNY2O_Ga2v2!kV4fCa1v^GfO!;9?>bxwpmgZbLt(wjBc?b^+bJ!u_LFPq{Nxg0C((?Ga4pD%?40r+IvYS9$hO_NXk>DO@o%iov`K%Y=&j`>)n)2C$vA<(P}9h?5=y{+KZKV`GGtc@S@eV>D~awG@~GwvXwK{Ms~li`(wI^qp;T$D^4-Z3Mu<Z-e`$QYq+%#xQf1K(M-CgYF2?MDIksx@gFsaIr<_cgZ!W-xSm>RXxLQpL^N1ZnyEd{!K`;rS5}&lvpF}zYAlLaBWF!$Z0S2LwVN1C0Yw)Yv<~Fcz>uB>l<rhG;3!>fNljmANvWaliDpO&=}Fv-bIqV%5NM>}W|z-!Aq1l2ix$cr=?F!BG$TKPqJ4{}_q|$3rVkPm0`A!-%-e#Ho-tTpnwOl3u}Yx47~GbPErorL9PuM{0vkF=BDz7goWTCAhUlMyRdd_VL?LLD36BEl(y{i8hpn^wSGF2rCmgwJI-Wm3YeD1?r=afCFkMs)uU1M@;x&&zYO7(bgM4x2$yLZHX@T^dRCssWWil8Gf33&kp>Ju~)`Q5OUreQeT9fS+Uey?3VUtLK(AW+;p8w6nNJ%JFmFk=_l)#F?OAW;Q;{0WrR6L~&>Oe~hryONMXgOc4B-W-X2rMwq@;Rg93taBTjBiOBz3N~Gio0lWGxMCn6v9K&=#hG!)3_p+)+V+VKN;@3yaTkw!#&090K2p68f~J{JGgZBa>~)aidmSL2#TI|85)th!JGEdUoG@p4>GQ>)UNIyAkDOX;J>VtA1x^W_mh7`_Y!F(PhKY{c1mg3+GwEhSG3(s4}q8m%}Mph+q77XlW3TNRub~?62>xTjzklQ0+f_faM48>NVze<%n0&+LI=ilpD+hxK3!Q&P&Fx2^OCfBEJ_4)v<x{=ai=K(q=2x8tBpxq5KBEAqRI12s2}Qfh(@36t{y`GWBd+U0Ifu`!&p$FHQJdgEN5m^9E5LE$)t;%c51-*Kz$-P8`Sv$udh#TqNgh>4==Mr<LOE^GKk4dYn-97<t&B{=nV`vA=tB`p`C8!ZC#v=rY!Al$|-0kruM-FD&(<3OHhDVQ1(v&s!2KU#+jjOX4C-rWp>^zUzKikJA@g4n6e^C!Axh*q>SO!$_mI7ytyeY7iFNUwIs62ViX1z>l`PwQ)mbvka7O2fq-bIq?Fxt6<}UoH&mK7INhW>)8Xrl`-;NBR^+cPU8F$Y6+GhPI$`OZv_IAU9<Exqd$tl-ClV6@ml`sTLg6pJTnTy}eh@cqZ%HYjG2sr{jeSeWyqdPVd3HE$vydWo50QLG>T4;7B5(L@Uqq<T{!sbc00komnVfwRxP!d~<eRdcd4Ete&ALq9RZx#alKT%-*Md4Vq$iB)6+K&ND$MTDErc2tlretf8_9`=<r<OOq3;N?so~x<s-8~j>(!ouQrwUuX2zB!D4KanieVTlHGa5KPiMw*l|X}nSA%{dc~H{K!vQO3*X%reGMUT*00CX+;aM&wW3DuMK40kiCk+N}$XIgDzEZQGhP9y+pk~?+H*YmFZsX8WRX*jUsxA9Lv&d4r9ezc3stwcXhI;O(i9a@zfa#?eWaA-%b6SZRRKcR_|6By6f1nr*efa1qDDH-)g5WQ&Rx<e`IX;4MXmz-Q>81@1@|6{c(zXbCf-D#~eSXWl-0YP5DQTgbk$^jG;9KW<yO+yS)>`%@CDFh*;ggu3_Qe9B&78dvoT2g<_jFafg+W*F?wRHmK`WWWU=q1k>nc>#*f$s&mzW74<KcnMX^Y9sDG=4L2XSztT|&7J`Aq(GoPORJ*|;qwM7xtjBUAMnrvr%3$A-@R1@gkQNBU;kM5gnq<Z(uFlUFIJCs86gCQl3p7k|iHqQ=ZVU0H!Lk@a?A(sDAVm+~&;gsW=kN(e;`^I9M$Z7+7&uvZU3ode4rYP<vmGhjIhhx4YYKg>;I%QxlZ&)lVK)CUkkK%*U34hqSuoP(y0V7~P-<z%ud48QI5O7hBP3<uj@7GnW4gs@neQ+W#N-eCh-NC!&gW~&__JyLb2pmd?@t{6qO1f|W;nuu&&GG1g!nke>d(pq3%6#zd3i;i7PikpdtLcAbTUAA^V3IUc?E}4C1B`_l7)k2cl<ZN*1kEL>>TdFBW)Vxp{9)Jwh0S;ufZNU`5B)!-@rl7Pu$$y|>Wm{1^9jhYzW&ju4U+^IKu<BNz-m6PP1bL?2JSZJPArG4EVqA(T*a$5_7tbvXgm`rMd%)t;mB2%zhnCSeGp#6>9J*S_{#QvY8j}Z4p<+C@g$DzxX5jqTQYhW^I_~f|O-ZHo>uf8v8{5L^0+JsF$GNsN!)O^j?sm+nY$?t0<M2UW4{sjMQzj;?PD$bH#Crp4_G%^i7e9GUKwJC_$i$TAPS5p{)Mb3;t3rUBXa$<xb5sT%>d1qY6nDpyd}p{Yi=-s}io^DCAD-96z_NGimIAXQCqTo0n{=9$xB|wweyuJMX`#koDH~Hz`Y}4$p*|-zJzMEC6mJB`^rGX3;X1f>Q@0rf3^@l3QcA$MFsX;316jfMDJLLgnrzcaZ8jD!TyW;MV8Dha9tC<8Pch9$Mw*g`oOotwQKLUN(-$!zL#2;8n0_{ltbh;_EJ?9W!n2hOr0vJ#ohip+FDOP}5<`u)PHs1=N~XC0meD&V0Pr^xG9P!8ZdZ7@;BH=tG0OjaQo#Uoi;97O3W8Cf;&%$_$f%>6z@$0ip>GEdXL{R*iESY6A#%1e52+A%UOafc{kMzgSf6xZEiFxp0|U$w3X|<209tcfCOy0e2W2PBf0*I?5RME^D;{69gv(_zwNU!4ZhU7B3k?9WT=f)=rUx2fpb(a3Heo0a(tUZCbf%3oTQFcXqyrie>9?kj_U$6qWnCstFn-<qM9xe>rfH<W6c1v?oF-MaUD|H`B4#O{s9Q*D6KokCi?hTlWx^n#Q6horF*ESgQ4fu$cwqaeP(b6T+e*d_&J(9m_VJ{X%PcCO1uUKZWsExk$dg4RADK=7z_-mFk^uFCY4;ebyrs+YF&T<BX()IafQsmrQ#ze;FdurS`(Wzu`lSU!#!LXjfC*4|?&R6RF&!xzbVkukN{65skUkK}DP3h?dAF5;jx;c-kN9XzfYXyjmF)mz%q03Xpk_d)0$>-Hb?MQ*a*xzjuo#|?6sQMNzI{l?t4UP`WlPt}+&EbmoKBFAj9JGu*`j44nqUrMExrMwVf+ql0Vy9_G76dQ%AICHMvJz6P9d~?Zx=x&k*K>>oc9Te)MnVtWRlMsTM@k@<l)G!-2AcK?1no1;}pzumIH;)yu2GGQ%>#yI4zxnW0nm56e!<QxL)w2Cn^q5eP9w{q#Z<>XU~`jCaa_#*~ok7RC0_#8b*4@-Y%-9h9{d2GH`<DIwOEQJS2`|OGktOBGZi|*^^D>aM5<8U!LL>PmD7Q+bWn0?4wnRPac2_kIEs{pT}t!GTmQ}WFY0fk>Mg}z&qS{3>0sIWe`kXoB9$L9@szyg!gt)#VugGLZ<vafx8nw$GnkVciscOEvYz+vI&`j>4vyqIQ|dm64<^$o9Iks@g`?#GSD1j_!#+q_XXzh5bl5?gj6qe^;-wf$RVR;<gs_mL-nW_`rnW)0#F4A=b<fc7LlSi<ds35X$^&<+k`S>-^|_&G{d%0k^()hkBty1yX7Wc_Kq^{@KlOb1R@i>LJUR=0heX(ASyOl@+n?5n7+ZU1D=DbKP5gL8U?spJZQ8`)<##0FP^oy9s-N|WD&01RLy5&Jbp0^7bHs0lE?rWMSl#VB@hoyRZLc_EA$YLHkE85tG>jGv**+!5{ZBfMP-v)d3a`-Stvx(bVS@3Ea6dvg{1dBm^6Wr0>^TvdF{HpxlHi&c8!rN5yd>}xpSx!pfh-2s|vzw#Yi|#O^<WwE_fVCG}0*=x)-juwC5yH`JB@i0Lnx4+yG{h+G_Q3F6CdC6F1UpaAOJ!r3gBz!lDnN@y{0~(I{8FZ;<WbG?I>dvZ!L7&3@MLTPz%aA{4JP{tGxTTpz7wKd>36eCkeu*y=+zhD(NkD;bw;OeC{LNu*;Mw^*{lZNWJ8hj0-JdD>@yB;0T<c-x0Z#cxALg5@wy5F#w^9C!n+ESKf&BA9>WKr>-?Opg~&R$e00s;S~GM;aE(=G>lM0@N(D$_?WpT6!2ds0wnccig@u9)$=WOW7fhhGFU!SR1ouoZ^8hBukH~(pzQ<EB)!@oR@EC(!yx6b~@3zH?K0tlwgNO$G3|rw==?%%DbN$vwGTc@}SthPbPC-iW}k~I%SY}2Kfxd)^QOpkPwE6xa{4lX*RgiFzXM%+#PgrqD&v+0c&Z8GbY-+!*i%;bc^l0cE$W?&;|9pR0`8Dpc`zMnKz3(#{+&!DsVI8Ar2-{L^-;jN0Y(_Qve8vXhj(l1Z{4Khe`<!KG>WU1ySOj*DfU@1=u~Cj2)NU_@*UXWDte@iAbvs)fQ$8SL*(&`=9PK50p!g8%<9om?Q!~oHWtfMa-R#=464HaVkuxm&nsW9I~qA$qszSzbXa}k*D(fM8HIk*awK9e(rIHOFDyxEgn2(+HK>YA*zP>_7o2)ci`OFn&oC*o@O>J#pBx%GxsbtiRh9|qW{oyA81#7vZlHcK=sR=`BOS$Nub7B354tccOY~a`NLc-57;#qBsCBTS1sIOx>_h%7Q)QY@XQ07tzBQs%qjqFbvQq@g=6qJG^VCDyfxGPc*N#dFgK$a6}8asA-%Wuj@HNGR?z%pQ3ZHME}Xd>$6{mwFAGU2;M1_gvlHOnYyIdQSu_u?p&`3m#B(p3wze-87&HCP<S%_Q!NA_eh-A<0J3~7dY$#tA*CTW^X>`!D6Yc6_w|GXbb;DKH=8^7B9-fk(EUIQ9AX%2y+)$mFLTwqM<*}y}9*)Jxap=0f^C+VMy%hZM^|T2MKTl!geNC?L<ocQNWhAW$A0=C60v%bGEgQauBWmzX`io45w+@0^%~#aQi5>1)hM7g8NwtwBpd8b<&2JXL9Z7zT-mq2|w*v%x#Hq)|SuolRs0nuz+Uj>$H2*p?`8}0~a3n~J?2s(b5|6>bmQL<0RyQ+ZiU#m{iCO_0=nIc*jUn+M`i8y(%xu~!bI<^7o#D~$63>)El}@+#WKG2(?{vx2Up<ARXOF3u9u2dE1AgyuPJ|WjpmFhlk>Kgpr%s_SEtvr77;Om$E5@XJBr!?63`89QUB&#TPGJr;ZND@PgI*(-nmir@skh1;^m5?X)2w0X;N~JV-wE#CE@BX9bAxlqYxXiXx6}1I%*HUqW4@S2mea!ejwv5A6NMA-NOU$(N_u~cq}WFj*u&&OgpADNt{T=HhqShbc!8l(aM;ikeS_#iPeB2uj<l`0EGa{d=Onbbd1NmE4MUFec9BH|F62Rf)6dMFiD2wuUif1;Hab8UPGm=`5W^be0{ZndcqGXfS4Cj5&t@_D$2^Zrlt>;b1|vLVDv#}JylEhLw1((W59Vb>C)1Cz)Jmi8uR{WW2QK;0UzzdC^aI8xi#&&N%;_oPBpl(wLPG}^;F539+=|iHdqzHzNn}eoeCZ-zHUj{Is+V|ho%`-Fa;d;yR`1N%1@%-O(m^{4F|><+H?ybqsYXNHe`6MBlpO>ymd+=J$``|@-Y%j!3m3+K*x9}dn211y<ncnWXDmlYA;D)dkPf6eo!6YsN0{I=fC3$+3oq%~a-<=8Wgsp61a!24I6OT~=)*>HPX=QMZ={upi-WNe%aNu(1)N8x5u6Z+u>mF!bpGneBH1TU1OtS0TAn8LVWtJqoEX>HL>Ix68JyC1NB9cidh+;__RyP^P%%Nf)b`mVi(%Cgpnf6GPg_VKGRgHK-+_lfO?7}4Ger{#Ha27D9^ANMDlMusQ$B)SGL8z6o-Fb>1)BR|d|Uax6pc{nhTLxX=yQNbLmIUq9|-Vd^ijB&1@aM64ls}}hlSF@OxE%qJYb;SJv@@>+oyB|GIt!+0eQ_01_yUI4eaken{admH=9Ib5NO_|E0Lm_&i~s*;PadD@+2%5F^t<OK|&oXeM<(G3ARkpqsvgzHAZ&>55pmPm~?8~%w_=LGR2GLFtK>YoD0qbIGjp`zto3(1`u>lxlX!9u${q|fO*|!Pfyw@2PAWLSw7lZ8RyLBJy}HAp8@94O{R+~rMo-=q%YM|K577vh=Q<sM|)bQbH;Of;zhn}!pA$@mP@{XCx(2mB^im0d`1I>M{>nUkOa(`s7413wl<g(jbmS4OC$9g_S{^{$C!H>Svz>!w~LG&gNOjahp8s%#ht0)=;2jUJ_yY+MUYYg@}Z55{OlA*AM$}!h;p~nDttu*v~(?))tPn=GMue@TE&SX9rS;B7^K_)N&t*8N;Xn)OUIvZfJk1<)7h=|crWS6B7jN6{>ih8r`?NA40pS>V@pSqX=%n1-`-W5J^3~eQ6BAK7Hla4Dc!O(XS7MNc!i`!^dz**s4QqxPVr_|4_Z8CpWHl--pF%KTHS|}!KaY$V&7~SsHf9?3T*euqH3CLz5l}vT9yvA3PK@3%D<#zBmod14<wd&gfd)E(-1%QE}O8D2}x-s{vaZmNn4t=O(c$Z*naZo5AnF76uU-p<C4ee9JskLws`6mPf+2Lj%oLrn<Pl=>1MrMWU4x7i=bwnVJvqw7-mp6qZia18SVqacyI4;iBTij@$dr(1M~sZ=pEO0i3ecXss#uqbsG#r;>^gfKH>$rQRF%4rPZh0W(OY`F<xQchZBq=LozAtOs6Fdg9ZU{Jy}%2wZYNvO2MMd7__XJd}ooO0ETZUtO0>UKt5VfGPKgs71b>slTYE<+Z6@<*)gv{z@mcnuNE)5AN>(9OkfW21pS*g&#UD|WZaH@IN{D^qFoPVdekKOK_f#+`eY580bac(;LFqO0pT}^)nVk<l&*u6tswAvMLM%Uo9>SH^$;(5h}<1klqFu9`7sjM)D&d65YMUp^hl5;#@*4hN6J68TcKu7dW{Pnw_ynQ!kyqJ9y4?6ekAMyckS&WVl|lYz^E;DZWwmy7pf#o7P$i>e&p^E(ZHAcuHUKH>&?*#)i1HEnf2-Et7Xf_tMmt>Z`y1ct?oC>lq4PkCU-fLB9?)>c}@!@JBsGJ3fm5o!T>_{ILJ{vb;<A#XuRJp;<|%rA<MCNJ;@e0QQC4ZbstY@3*dI-oaDEd9|B`R`*5fa?E(}Cx&a;DE#K&~pNz;J4-v%=6gfD&J>}y;;fw3|^ZiDjrYBmPHR+U((2YT*H0CeZ2dx|+pk%S$E<(!`j33u8tpl2U{Xm1I#ao&@`Jzug>xz8vEO_AfkUae%AAyFvYT}?JAU?x9p@2B)KzFYrWqYa0M^FS<#^Yp+_q-L4C8y7dO5%=bnWYU$G+*0{_%%~CT~!P*dAo?bR<kW|kU0gU3E^NS2z>)nK(J|^e!{R#5)y7mgPU&<gDq}34;>>OV}QYOtm&kfT9_1{tIfqbNOSUW6_Czdx?Le()orTI8$#`AgxN|anaai8ln&b4eGqDq^689he6t7uG39+;-FZTZ3Ian6xpZ5=fqI5GH6{y4r8Rcb(bxPC5DGH}4V2Zg7dhPX1~LU0fN}G@QezPNc?0#LJOBaP-V`}dRA=-ibkDhL`MOExC`Fby!Rch~Is6$^-P=X8?H=S@S~*)!FH7+8)4sbjC~Y4ABAM*lEgsIXNmJx)_lJ0l+E&ocK$9gN8O|vmxHNoH#DgtAJY=K&6wE~TNOlPDJs#bY=+0~za^qWaaou1r>9C06wV7`hRZqaM-^W}W3r7PQfJeYcJD~(rB_`AT9dqws&|%~W9mC-i4u(gC4!FC^xZ5;M^+L%OI!{Z-v<_x@odR~IUeRh?d3K^gh#E1{?X1a?dSuo|@fblj#8Aa_CnppA-YoJ6TIiBG;FD?GX+LGCBQ4#lJd}tThWYTI(R*k5%I5qi(fs3V0u9merEzbW2UmNoUAik%OuzD=*@pKX;zK^dBmhet-koAfPaI>Gorzo+P~m2zS7$Z+fL>($lg!22MNH_6gE*IFg{eJBTt#yZK;*mRqmfln7ryI2dzwM25z;~zd1_A*CYi#`sC)kkrC`g~DG$(3$BbFp&&$O-i+Ci_zGoxvTH@@bg|?uCorovYra}lq?pWQG!A!c8Zk{Xx;(*qscde|9Al8EnCrmtHA5d(QIojv}lYOVWXFNXHpojUxmTULv(A^uj8J@{WdgPT(ALmNnF&@X0Rw*VL77A5G+X(We7(Kf+%Sg=isxEo-+4d4>s$^Dgd9udH8SRwOJ5w@JoAi^=$+wRt5d?&B*ibw(>H5K-CuzB}XB1n;{V5Fg+iC^txv1SEVG(tG0H=m<$Nms48^9Pr2{JNL8Y-D+n%s0d4LdzgaMpNqlZ2xO6Tml_j<;*<Xwlrnu!QASlZiLmDMgNbGGTIL?uh+P1x36a>VFk-<cE0B`3@#Sj%^e?vd`ib7m1~S-S5uK@8fbXDIJ@VBBR-F1RJ;@Bmqi3Glul4dKhWW6M~#!d~5DWu1^+udQoBqZzM<Pj~Pq_=*IaIE-V4zfHxC7RgQO=;wlQao5}MbU?h<U1^NcFuwY{B;I4&9CAiTYb_E8#odQ<(&Oq2D4=~(lMJKaPjX3Q39{AA=4+xcsFda{W5+HJtKyMdSy6+gi%-p02v^s|SLpqaUsd+5{QIrtH(rHUJKAI<^!Y6eSy>z|_w-$$f9Z)8U&sqWZy_-gw4tE%Y$R9lfG@}-&aul+j6(AHIcRkWZZxV`Chy*~d&Cr1aDWEzUL!T_FX1oMQobx6&r7;g3NZ;j}i!kSX(D4k{NO8G=2%Q`#FqOQkd^D@D^LPCFB_V^^Eg+fI3@0`MGC>zUk%xe;Mlsk0-|z<cIgQ4&J+FZi>>&7)3Viyo?qrzhFkXxvvnPvqzCys5J%4!$NP=ZF$2n_6y<1ze;p`kOc)6MZX#lMA$r;fZD4{*|g<G*h#$%l-=791LQV9kfb5*adHKS4~Pf^B_*}s8pV!T)cLY#N(kX&{lyOrVs$M-oNkOKOXMb&N8=si$XQ$G7D8P72F?O790GI7aOrtAp#@Du@A&hZcr3d=zKeYJEjSTW6-lf8WtMiyFrfhyI9-&x2eWKV5U@OdNjryfm;{4{;d)7q%#R;w3iy@2MESADWVq6NB_b?u?onzHF(<o}u!qJ*tKiT4O2C^XqC4&c!|a=kxpUo$wgLdn1K)!EawEnGLyj`n><M*|<K>2Wxj>qGhk&v^73lh!&|NLi&AxBhhgmKJ|dHgoSP5C-S?#*;+`^SOl+y7Z+#VfG?^_wC#AneI;?9Sk7A(*r9FncZdcfNsl11HxO5OExW&<wGlj!M1!1Z-dMClrJ-rs?+hO)2H3j0}l6u^0j9@1qw4)X=hqVZxRky7qvlPzF7nsm>L`oLmg`aDUx)1c&+Tr3)HLTDdTtUbqSv0Ne}cMSBii@pc?##tDTj>_H&GYqH7wU$dABK2$lSJ9L~(tQ|cPNTsPzAva7jyxxa`9p$Y+hWNMwwpUW}?m)zGkiz<Z=a~gS0p+Bb6sK!W{2v26-;&s#eVWZs44Y?!J0O&qw5AmRUB%O6E)+nN}I=SE*#IZobst1hYVlt8X6i{mUnDe=z@%C6GUO_(#JvifxK*1=1kg|ZV6&MN1B*eFi@IgS5HM7KQ&ap~Y2rXEB&YnU*2n$BWqcUk4&MBN=Gf;g92*-i!d)^E!B)0M%crgq;pQG0qzziVq6cB_K`8{(l-f$B<BDgT=ZyN{5lbqFVp43cQqPi0t8jb#BQDvF|pz=-?a*9VuQNjVmS)|v9S52u8Ae*)wawM`1ESE7l#Umq$ZHOE4n0*+gbSBuOk+aZxO%~vB`)HsEP-b{|^NkiCZ!rx&ap$b5J~i@YVOFBs9U3$6&$Nv_Sp+DA_*8E4In53hK)4Ah$Zaw{Ao7q(njX{KYJ^RLvbXL-eTYY;F;M{D>a&<2(DK3OgmVm|OE9P*Kag{@7Z4crIg#aNzE1EQjHPZ{y}+wlrw>)y3fu?0S<qBHTLcvs+z^Afc9GslW0+EcG|)zhLrUb)T@)T^Fr#@x-(sgGUj&58q2ik?<_^*p^9;Z@MvFxA;CuQ7s9K%!MKb<|`9alq143n35!i{mH;_Jr0Z>j7PF9_49$nxA_1PjCp`ZuoxMIu4>yuw7x)zaEcm&3#0>8BX4-YQxCzKTsJwN2r+tXVBQnke6y{s`T3BW=mJ?Y_X#3PUW11e0}2&(HW^+LeNNETX7VmRd^D@Jy{xv0{q98%Sw3B5mA#ApeM^VBAf#lzWY7^u!@&>gkq%O1a;oMwQwS*L^u+QjCkd{BiQjEbroPR%f(siPWnv2Zq`UJkc+AMIfv@Z!j!dgi@MHI2y>^eST*cin-`VNf=8^+rKm;W0yyK|fn$v~@J|!2(MhFd(#U&ggM$Jb{T+n2R2>$5s$m&^rVVYRC`ynkQ{A&ujUfoYXbJVu0%&NL!3FHplBjJ|-T6S(ua8;8<<Zgk48^gtu-dAo-ui;_wQ1#UROX8s97e>ovef;GnIKq3j9jBT<mq-X}!}#^N!n)jaA@OL^P_I^`jtUvIwys9yq7`JEGVXd-}4O3tep{m3+&^0}igOhmWNR(2l_F+RMUr-A|Kt1mT+(v*)<L@91Vd$Oo{&VYvj0JWO(=~WkZ^wXQM#RH2$cByzy1w^c9n_*(Qhj`KSiLfp0{i}=PyaJE(ur<VE`}pgJc-2j4Gyu{4KplvZ0ZG*(@2VeA4ma!3MriSzc<fMd==8l^1bI)w<uiP7@n_M67R;psao3U$taYTPpw6KW#BP(*qh9<M$mevu@3`+vJ{`#=AW9s>72O^*B7jc3?4OOaGa%%il6lCX?9Q9@-3q3FlBCTp&HV#}Vbp4ufq1*fZH(h+tD9yDNd7CJV$B4(1q3VzV02DmF~=MB3UBJ6KnY&+;Xyge!WUBh{tzZtfb%3ahaH^mhj!UP%(g?BPS>*zIILOKLm<;BpuDp|b8t^(vIQnAyfJ=|pS)RQY$8}qMlCPoOl!v01IRkg!ZbXBC?2|J7+XavK7|(SaeT<vXci`HITh%@?X`Txyxz1VfS4f|d35kamqsFvBb$6no*o+}N7W-ff+~!PGCN}Zi7P`eMC9~j5#Pu-CT7(xju<*#nulS6QeFhKR%*~sXNFZXOiqS^^EoFYKIAJ-Slt=jH|4`>&VxZ;f{RoV4_VDS^N_C^gE?cHcJpvCb86cQUh-9v0>PzCcS})UGg_V^%CkjyPiQ<iWp?(3^fw*mu%P$24IrTb*`EZ1y(}Ds5p-R~?my%s(?hv2tWiV}3Z5mN$Bp^=RnL_w`Qd;%xNAmx^POOl*saRYoZT(ez{*9}rD~W|42FAHuk`R_5ekH{18}Oh<;TKqhb1t>!fD@|+vJmHmO!fkwye=<bX=vQp)J7vSoiJ$;h5?rp#_{)xkD1l&R%wJ_T+SmB$BS-1||-9&bzgNUq&FC5aPF)I?K3bMi15}i||uCYQX8-#oz{39#hARvD5n`2S`o9T_%j|q?v~@C1aOP@!Zr?hMdA(6LWQUDi@MR(5_~hb^s?mLL=Ry;m$f#k9jCTe#7><>$TK^&~4<CXS)>hX^p;W?Dm^ACJyG2?VT;P0`)ZG;XoUVQ?wcO5eB*33J8lDz-2yj2MZ<p!sywt!uD!ps3n8(3{{)VGe?hYP2nLM7{LZ|t5RfS?lXKsIF_5woS2xvQ9)WTIg$G|lF2hGN!z!J=!xx~j0LsSx_7j3U{50X-<4RR+SlEbMqJt;`#elFkbma4BxdRa+?wsT(8e;;x^o%=$B{lj5Pb^A9`aGV$IZ!}Fi3|4>@`q^#Eo|W<j60IizJ5gP#(!99o6o<T~s~P(}M@#_w~oXXH@*MXx=~Nt90^#6qLb)Es?SmqufD9KjZ`6TgqKlXUkX9GJA(-`!N|e!}^C8`|zs+T%4uSl^^IJr07C)fj#PLCO6gWN*5~GM>o73BzOQgLCz7LEMm9aJw%NoZmzxL@i=fkv|SO)bJ%0Bk9zxu5C1U5p$=I|?`$(p^<d_33%~-ZgvmU4gT%w-IMO3x-o`?rDizcBH?yJaG?`BtK(INk(7;DPw%xI9eX_`e*^DT1@8P0$BHs&zG5tDWOC$e>Nk<f79?}@z*G+T2r0Ho9mY2!?K`|iT&7%p^#wiYQ2xDNFz{9v~rfE@i1b?yIB{UGI5MMyB0)0TJ=C6pB=0C2f$An)n@FXrqPm03%m^0<__L|{+Gxv}K*)5Jflu#9l*}XY#l(IND=rHmj4KFRkN0K3!orXjy)>FrZ*X*7lH23alz;`N`aZK5bV1o&0%}{e#RNgFdhc|v1mDaV#(0@SW4^n7jb3nh4f_0AkAq`0E0iE@3E@>Y6Yiazg`X;l~lwhRuLmEaPO&)5<!+_RPR)?qZSoN%75F%~~y<rlj-3>1VG!6*m?IOHf55UB)Thmb?3Bt9e!!Vr_Rh*g!&jnA^DUEZR*&Q28nlRrigG=F(?B;WWCY)I1Aq=B?=!1bTPlpD7Pq&vTJay+sfJiXe9ds=2y`+H-j2D{W*?7CCvQNp<F_?Es12q;Sj<G6h#G(1V=w|`8{ZIpr#1)Roo~>*oVc<-+76W)mFPpawB*ay1H#d0DUl!1O$pk*|>y5BpAr(gagXqD4=z$I+7>InFTsHm9A`SzsL9iW`<#rhpK;D|GvNs2}&O=+l@txw}0nsROHj@KXENg}y456U1OudrLzy3ID9K<#W(P3yM9UmzbP}$ynyL48P*dfsO&1r<`ATLPV^(TuebnxYv!y#7hcBMmrq0xlc<=h?|EeNv7JDl<msAhPJHU!wil%^uU2c)EhJP<>63a<}&%3;7FDpc_R9W|u(tHDV7xWPC8bpx|PWOaJIQ^x=eiUas|4HRL_Duj?RuX-Ki7AOO)mTq6!V)6mz4}s6)6bI~{9`Y`ey!09<h?zYtW_yA}^~Mko;_pjoEi?}&=XQ|zxRLARng9{0VQ&vNw!zE*!IH!w3}@yX+8-z3%_3?SGv>zf&X2K+lp)z6l;4^Ss&VKK;a_>k)47!(|MD_BbkfY-NoLIeU7mOfmvJuQFmRg5P6u}XhUF!Wg%q}|TTBOo%Yw!brYR1|bdQ&Z%jC@>5EBssQ}<=~Q;?(LlNrdpgaNq8JNSQX-Aj_&Fbo9HwJ$RzS^xWAEO<>)a>ZU`k<>6Am?Vf#qx<|gYIBrr=05<NZc7H&nlVZ$QAZS45+L0lxH-b1**L_UVUoNZid^E8r3VC8rqb(|kvttaX}AV;_k9aL7SV;h;F-u9_DRA`sZoc_B^i*kkh755+)YB>=6Ig!ehPuRIcFYg=13D(yblNfpqISxk2#(v^Tb~7@_jQDK2A5}4#OU)sY{b-mDr=->emXa^m~ypu0t&$8`RC=PY{yp2#o+Izyl%%I&il+PQ$SXJ9_kKsJ}hMin@pf7dK~>qzQN*bEtS8T9ERqmtvO&?jE(`@ro{_&JHoKVqDvvMGJ3<B<}a3-M|!Bg9*)^`BnU*b}c-s5MmqlRD6Xu96`-(j<(5hl9TO@;Jiodw|NOW+Z=NyTm%)(JmwfMKv|Cvx4vkRZ5<tl<UIFhcQpCc3s9e)REN;gtE{7a_*zsUd_{9eT9$3_hivJBg%0C=40Nt(G_4N25yxXsow9mf+>~B)I3q_fMFvxJTIWE`A{w99la|jBkPc)<=KESQS!hla{nTyF8oWU5IT(UxhXi9fIGX6+i*^VGngXb;I<ICL?XtTypwMm&1*H0KHI;8`WKQ*I2AhX9#0k1NbF^~~eogp_bo)7Mic}O>4m1+`3uHM6!Bq`p+?_eR4?OKWwV>cMCQA7n<LCDxgA*95O$3g?O|hw41LF7|u3$^vePp2FF(<PO<l9Hdnbc29-NKNPxx-YHBwzb8q-RVTn>J=2?WhW-8Z${bl6v=tIbPvCnn`Lxh$o=$QBMiAz89hKIx?bSY@}~beW{sh0mYh*VGs54a(Jh1*=^4rS=c{gW!;{Y$I2Qg*yJT;k2DZf!*kjy27-(|$+pp<R(MOKtA;y|)8+<m+<iT!?X(%>BQrNqLHk-1hlP!Sib)T0_y)>jQsX}FXS5jy#?j=qAA3YH`=voVoB<N{_E^Q8nr>u|!T~Ot#~w<|q6y{E+8+c&dZpkx*O<Q#^H$1D3o#7)W4Ig+y`F*}i)4s$Avle}fVuH&d&Ku3Rfav8AJldENUfu@y-a|E>#*nQ<{%@ZKh7EoEJpYrYj}4PIpE*xOT)=<-0D!=A?SKEc13WsJ+3XjXK!jC-)j^))!fUiUPBHKy$MvgM$!9D-e^GCwE5Gsi8w2S_>+ovV}e0#)&O5_5RU4}n-a{2HB+POdg}yjuRsBnXx;(!<!9j7jnk$0VfgwaO$NthuzERuEHV>F4+w(1L%lFTZrZYx_{N+M3_ne>C3NGrHNll%e-$nbSD+2qz5!Bzs2{ijX^)Y$#}yC)K^aEONdJyrj+gG8)<iQW>TN)G8iLoW3|Dl9JP?NbS`>meU$Q1q^N-?n48>anF}pP=)p_q3=Kr=vO?No~V|t9EV$|sv_t{my;-||?20r#|P8zT?zCb7PmD;R=NsF_)HCH4pF$<eRs~7@H(Yc+CxA@1Rz}i;<8t6pdW0_FJFA|w`U}1H{B(|sIW7tFUpnMuu&ykFl-+TA(V8a9@e=d7~P>`-a=4`rLJ^@p~OG!dy{7RRln`78w28}rOR0{ST_nkQe%1ir?MM_w2dMmj3di7-W>ZBY<JZW4VEqq?N?!fyo$5`yPP#o#WC@r?O1cIe;3NAX`mA+Rh_?RQ8#jmRuSzmO^_NUZ<I|J`ZN#@de#MD^HPCalmCgI4jk^NY-r)FntOU~23l{GpKIB1Sh5RGp77QWl()zPqUPa|&Et%>YL)`S7yZpECiA-z8i8*qmigd=(v3d3i8BCJpcysP9Y11sJFCcVlBWioH!F3Yq1ShOcIx`IR>jWO^K1^azRlKdX_TuF|V9gfGI&8+nP=l$sRKr*!L(O+og-Nmy6KR<tL&z#vOiN3i)uw>DT6L58y=NifPEyZ{4zTEF$w8r6(T#fHVxU08N|8DCxcKID9DOy4aE}A`Jk+66x2|o5rs@TRayEz*8%@<%wksrx-ba;&xKGq~2R~*8l|I$!e!h#asseXqekT=oQBKKGK?GM4itQA0k_<NDz$&@uKlRt)M)3e>+$IO28m`!~&1m7m79%~|qMFoANZ&;N3j*=myM#{AzH#Pt8Y0q?qqY0LbJ<}CxO@yl*pPx}uoRs10sNE?iY6eJbECcwl2mt`VQ*GdMPr%p1*mT<TcE)Kl%Jz;>+xl_Zpp#7_PGj71I7oz45T!EbzbJCibhxeY8SkL&1xc@_2rpY&&)!_qM#I@2xFC#e8Stm=?WrTt@MDqBQ<Fj}a*XC?uz3Ko$hyW^Q-GL;uE}GKR%1P*1Hy=6WCY<`nk((BsfQcKdCaNMPW2MW4#Ca{e-cPo@?th~xMSeli8w}?vl|5<H#21<@_Uil8r>NUwHs)S`m*I^>oP!MTc$cqo_u<O074YO2S<WHX&bPJC@w5#Gsn!%M@_dJ9rE%mI;nOi0_173ubp%(O%3Z91=zc&P)7u|UWo1bv1rew+Cy$M-J7#)NfBXH9%GhBiwMIXYmf-&+>VPQIq79N)ku?ZF$MzQ{Si#_(IBfcZgHokaSb{U>v~@t?a$Q|uWq)BgD~JAbEUfSy+}EYTEyH$3WhboIG{Ltj}lt~@(AD^AA3U@_6&+{m_s3;GpGhg7BtR<=5zksnY5oHh4yAk#qNlOW`oO9{(FD^|M%&yJo1^pW7fHr`zN1YC_}0h31LqB2inr|-T'}

SINTON_SCHEMATIC_PNG_BASE64 = 'iVBORw0KGgoAAAANSUhEUgAAAnkAAAEmCAYAAAAeOYYeAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAAFxEAABcRAcom8z8AAI41SURBVHhe7Z0FmJRV+4f3/322n58BqGBjdyDY2F0ogond2GJii92B3frZSneXIIjSsHRJg0HDzjz/6z67Z305vNvL7szs776u59rd2XfemjNzfvPUycoSVQoTQgghRFoRzuVCxBIOHCGEEEKkNuFcLkQs4cARQgghRGoTzuVCxBIOHCGEEEKkNuFcLkQs4cARQgghRGoTzuVCxBIOHCGEEEKkNuFcLkQs4cARQgghRGoTzuVCxBIOHCGEEEKkNuFcLkQs4cARQgghRGoTzuVCxBIOHCGEEEKkNuFcLkQs4cARQgghRGoTzuVCxBIOHCFE+rB06VIbPXq0tWvXzn744Qdr3bp1Shjn0qZNGxs6dKg7RyFE+RLO5ULEEg4cIUT68PPPP9vll19uNWrUsM0228yqV6+eErb55pvbFltsYWeffbYNGzbMkslkeOpCiDIQzuVCxBIOHCFE+vDOO+9YtWrV+MBPSatVq5bz6q1atSo8dSFEGQjnciFiCQeOECJ9eO6552y99dZbS1yliuFh/Pzzz23lypXhqQshykA4lwsRSzhwhBDpw6uvvurCtKG4ShWrWbOmffXVV/LkCVHOhHO5ELGEA0cIkT5I5AlRNQnnciFiCQeOECJ9kMgTomoSzuVCxBIOHCFE+iCRJ0TVJJzLhYglHDhCiPRBIk+Iqkk4lwsRSzhwhBDpg0SeEFWTcC4XIpZw4JSFTG14mqnXJdIfiTwhqibhXC5ELOHAKS7Lly+32bNn2/z588N/rcHff/9tM2bMsBUrVoT/Sju4lt9++82GDx9uc+bMscWLF9uSJUvKXQQmEglbtmyZ/fHHH+4+lxSWkeK59CYr73MTqYVEnhBVk3AuFyKWcOAUFwTEJ598YjfddJP9+OOPsWKkd+/e9uCDD9p7773nREc6wxqczZo1s9NOO81uvvlmt9IA9s0335T72pyIx+7du9uTTz5pPXv2DP9dJKxj+vTTTztBqsk1s5HIE6JqEs7lQsQSDpySgMjZYYcdrF69eta3b1/ngfLgvWPdyp133tlthxcsnbnvvvtsq622smOOOcbuuOMOZ8cee6xde+21tmDBgnDzMjFv3jx76qmnbKeddrLnn38+/HeRcK677babff3117HiW2QOEnlCVE3CuVyIWMKBUxII1T788MO2zTbb2I033miTJk3Kf/zee++17bff3m655RbnmYqGDfmdcCShzrgPf8Si396HLrGoiAzBs4hHje3CY3kLt40eOwxrRo+P6MKDd/DBB7sF4bm+d99914444gi75ppr1hJ5fv/hMUIIYXMPQiHG8fDE1a5d21566aX8+xVeW0E88MADtsceezgvY7hvf0z2lZOTk/+4v0fRe8xz2W716tX5j/EcHouG36P3yv/OdnGvhyhfJPKEqJqEc7kQsYQDp6SMGjXKLrzwQifo8NjxYf7FF184L9RFF12UL/w8/H/w4MFucnr88cftu+++c7l9USZPnuzE1MKFC93PF154wYmdXr16OdEQMnPmTPvss8+sRYsWbjv2HxUWiLKRI0famDFj7Pfff3cLpj/22GPuPPE4euLECNt/8MEHziO55557uvOdPn26W4/zxBNPtOuvv34NkTdr1iy3X0TaE0884Sa4adOmrbFPjsP5sN9HH33U3bcRI0bkCyxEHtd8wAEH2Isvvui2feWVV5xXr1+/fkXmNz700EO27777unP1Io9981q9//777r6z/44dOzoB7lm0aJENHDjQvR4TJ050QvaZZ55x4V+EIcKtT58+7h6//vrrbn9ReN2495w/58kxXnvtNbfPqFCEuHstSo5EnhBVk3AuFyKWcOCUFCb+Hj162H777WfHH3+8CxUed9xx7vdu3bqtsS1i66233rJTTjnF2TnnnONCntddd50Tc96z1KpVKzv//POdd5Bw6BlnnGFHHXWUE1UII0SEp2vXrnb55Ze7/5911ll28skn26mnnupEiBcwFEogfC644AK3zyZNmlj9+vXt0EMPdfv/9ddf3Xah8OBv/sf+N954YzeZ8lzOD7F05pln2g033ODEEXTq1Mn9n3M477zznPcPbx/HJD8OEEtsd9lll7nraty4sduOe/Hpp5+6fZG/iDjCG3f44Ye7vEf+T1j8hBNOcKKQ/RSEF3nff/+98yr++eefTpSyj9NPP90aNWrk7jv7xhM7ZcoU97zx48e71wLheuutt7rtuU9si2B9++237dJLL3XXfeCBB9rFF1/sRJ+fwBGD3GPuKfcM8c+18Vp//PHHa4hh7i1itSjBKgpHIk+Iqkk4lwsRSzhwSgMf4IgvwpnsEo8XAiPqdWNSRxgddthhTgAyOZEzRhEDguHqq6/O9/q1bNnS1ltvPecdRESxLUUPhC8Rcx06dHDb4VFDtGy77bZOSFHg8dxzzzmBs9dee7lwJeeAdwnhybmddNJJ9uyzzzovEwIE8YZX76+//so/Vw/njEcLb9qOO+7ovJM8b8CAAe5YiBcEGKIMIYUAqlatmt11113Om0dhyrnnnmtbb7218+zhVaOAAwG4//772yOPPOK8aXjW+Bvhh3cMAffGG2+4MDjXdttttznRRxEL+0KcDho0aI1waxQv8iiIYRsKYDjXffbZx3kOeRxvHPeSY3CegBDde++93T1B4HGfKP7gXiIkyEe888473euBUKxVq5ZdccUV+aIbzx+v23//+1/3eiJa2Qfni6BkTERBXHL9PM4XhSFDhrgxwL0MBbeIRyJPiKpJOJcLEUs4cEoLYUoEDbtEUHjvkAcvDh4uBB0eH4QBj7EdwgHBgChjcsfbt8UWWziBl52d7bxRhE3vv/9+23XXXZ2QQzAhLBGCCDzEHBMJ2/7vf/+zOnXqOBE3YcIEGzt2rBOWiLT27ds77xHbISwQpAgW72nzRPPMODbikOtC0CGc8EwhehCpeCgRg/fcc4/zjPE3hSZcH9shkq688koXGub4CDDuFV5I7gNildAqIomQMoITEYag8+FgzpfzJgeQogq2jYZaozRv3tx5VvE4cq6Ep5s2beru79y5c93+aQGDkEREIgoBAUqIGLHOczge2yJaN9xwQ2vYsKG7Ls4F4YznlPvMPQb2t/nmm7tr85W9vE6cB68bXsKoB3L06NFO/OLxwwPKOOC6EaF4dqdOnerOF+8m+X2F5WRWVSTyhKiahHO5ELGEA6c08AGOgDjooIPcBzsTOmLGe5oQTHhp+D/eOLw7hBwRcwg1xMN//vMf9zi5W/wP8Ybnz8PjH330kRNMeN4QlVS4ss8PP/wwcja5YWFChoQ7OS5CEY/fkUceuUaOILluhCPxrPXv33+NfURBEPkwsBcpeKmOPvpoJzARYYhMRBz5d3jj2rZt64QYArJ69epOxJCzxrnwO55Bnk+BCtvyXN8Xj/0hZHfZZZe1qmvx/iHEEEN4vOKIhmsRatwPzguxyv349ttvnUeScDFeNzyJgHDjfhJq5Tw95OXhBSQU7+F8uQ48kOPGjXOPvfzyy06AIvaibWXYllAz95DcQw/3knuCuERwf/nll+66EZUIfF5DXmPOlWvhvqZ7lXZ5I5EnRNUknMuFiCUcOCUFAUdiPV4tRBDhSyZ6PGc//fST2wYPDGE5xMl2223nRARGeJNwX4MGDZzQInzLtoRrKXRgcvAgCBBWiAraiyDW2Ad5eq1bt46cUS4IQc6Dc8CTRzsXtsU75EFckONGXhwh2IJAiHA9PB/PElDogUBEjPgegBQbkMuGxw9vFl4+PIBMdHjgOA+ujzw/hAyhawQxguySSy5xQhmvFQUnCBu8jEziURB5hD8pxCiOyGNyRei1adPGhYMRu1wv9w6RSXgZcQ14zxB5hIfxLgJijZAt27Kd96Zx/9kf23tBiMjjbzx30Ukd0YpXlfuFmCsoFMu1I7zxcOLZJIxMbiUFNYSq8eTipeT1Jz8QjzDCryCPZlWgrCLv3//+t2200UbuSxaCPzT2Tfh+/fXXt//7v/9b6/lFmUSeEOuGcC4XIpZw4JQUBBCCAfGGOMNjRAsPJgRClL5ylrww8r0oRGAC54Mf0cZPxA2PkfgPCBhCqwgpLwiYyPHYIfII8XEchJPfXxREzd133+2EIgIUcUVoEaEW9eT98ssvLlePcCHbFQQhZS/yfP6ZF3m0iEHEIE4Qq3ggqSom9w1vphc+5KixjW9hgpDr0qWLEy08D6+jD5Pyf+4BXkueHwUBR4iUyb0wkUe4FmHH5IrYI//O5z6SY8i54TXjHiGIwYs8rsmLYUQenj6ej0fPizzC04h0ztm/bngXeX3wFEaraTlP7h+CGk9iQSIvDl73YcOGufvC/hlbt99+uzM8i4wF8iO5xs6dO7vXkXOrKsKvOCLPCzlEW40aNdx7ldd99913dx5aXvO6deu6vMnQ8MDyGvO68qWD5/F8UglIqdhkk03c/sNjSuQJsW4J53IhYgkHTkkgdIZXZcstt3QeOcQUMOkjgJgIEBQIG8QggoyQXVRo8T8mZkKzPreL5yCWEIFeVDBpE9rFQ4UoweuDZwlvHaImWqWJR43jIIYQVuR+4UmiSjR6bAQH50lI1Xsd4+DceS5CxYs8H65FbBAG5VyY0Agx+rYl5K6xHZMpXi/uC6FNvFAIKv4PeAIRK0zWeKt4vq+uRdhEQRQy6RbmycPrhchDyBFGRthRYEGFrZ9see0QkLxGbA8Uc+BtxQPpcyoRedxvJnxe66jIo7oW4ehFHoKD5tgIMdraAIKOe4tgRRR6T2hJYT+E/3md+eLAuVLcggAlZH7VVVe5nD/C34SL8fLhqeU88UryujGGCipWSVe453jhQlFHDiUCDFFH2J/x4N9/eJkZ83wZ4TXEi8xrwxgNjcf5P8b2PI/n47Vn/CMQ8UaTksDxOG5U9JHzKZEnRPkTzuVCxBIOnOLiq2WZ1BFTtNKI5uAhMPj2zyRAKxUS+ClMwAvApOwFHct3MWEccsgh+SFTRB7bUbQQFXkcD28CHhyORQ4ZIgPRQ14bx2BCx7NIIQehTbxseOwImzIphSIPDxXtWooK1xKORuh5kcK5MWlS0ICgQXRxTCZFX6nL+TGJIt6YHNkPvf78ZIvnCfCaIRYRYnjXEFYIMPIXCdtGQQRyXwvLyWMbRBVCh8kV8YMQR1z6XDmKG7jn//rXv1yxCCCcENG8PlGRx/8JERO29a8Hrx8eULb34VpEOIIDUcH9QayyT+4BghVPb3kVT3BdXD8Cm9eULxjkVfp8Q8LhXDcCl7w+X9BBTiCvT7RxczqDyNt0003zRRXeOu4/70m+lODBpq2NF3SMz4IEXXHN7wPhxz4Z26RD4KllTOGF98KTLz6kYWSauBaisgnnciFiCQdOcWCCxBPFh7rPnfNeKQ8eF7xT/J+JAIHDREy+F54tPGh4BBAPTAz0YWPCBiYuRAFVsn4iJoRJjhaTl68GRURQ3MC+EHr0b0PMIYIQX0z+nCueQlqCICajVb/kxpE/x4SFGCkIRBjPpVDBe/I4N8K3iDOEJN4swtNcGxMrRQmIDPaNoGMSxMOIIETI+SpWQrWcM2ExRDCCCRHCPeA6yEmLglcNkYi3Kq7tCxAq5j4hfJlcEZN4MjkmuY+IYKp2EWmIZtrT4E1E9OKxQxT5Bs54AvHgIXApAvGvB/eW6yWv0At2wulM6ghMrplr41x5fRFePs9vXUIYHM8tYprr53VC3DIW8YJyrXj7+ALAvSWfk4rm8l5/uKLg2vCW8WWLL1SMScY6eZc+55XXCeN1R6R5cVZa4/nsh/35fXMcHuO4jHUEJl9SCAeTjiGEKF/CuVyIWMKBUxwQTnjHmPzx3hTU0BZhR9iRXC4/wSOYEDmIICYEhBATcbSnHl41RA5eFy8qEJG05aDyluR9D5Mz+VjsB48CExteGy88AMFC5SZ5dL5IAqjQJZxKvld05YsQnsNzCXd6McC5ITrp2eevHzGJgOAcyBfkPCgiIPeOymBfWYpXkucitHwRBF4yf84cgxAnnq9QfCJe8GjyeEH3nW2Y/BEvvFZsR+GLbwTNveL4iDp+InTwiiFA2TfNmr2XEI8Z6xLjmcNT5l8PwuUUWJAPh7AChCd5l4hVRJQXGVRR+4bRlQHXwtjhOt988013r/FOElbmPBmP3GuuEU8grxlhfn9dqQjhds6RL0e+6AkBhueS0DWvM/ffCzPvyStvi+6X4zGWOT7nwZcdCn04Rzz2vA8zwXsqRCoQzuVCxBIOnLJQ0g/wuO3jHoPw8fBviMv7idsO4h6PeyyOorYrSWiKfYXbF7T/8PHw7zjCbbhHcSHTcDtP+Hj4dxTf9gXhV5BnrLDnVyR8qcDjh1DlvPEwk0+JOKdfH55HPNR4XxHn5FLyRQWxGi7RVpHgvSXPkPPGM4mYw7PGz/Lw0pWXRc+JcDH3lS98fIGjJVHce1UIUXzCuVyIWMKBI0RpofiEsCHesrg1hlMNxC5pAITbSRXAy0wlL2FeBAkV2t4rRjEKnms8wuRa0nsQDyliJRTp6wLOlVQD8tu8t5j+kt57Hc2RSxXz54S3mrxXzhchjYebdI/CluYTQhROOJcLEUs4cIQoLeReMbkT/k7nFiaEQsk5RIhQNEQ6AN5J0g4oaqEXIsbvFMjgVcMruK5EC4KSQibyMQl/k1dJ3h3CibzKUFyloiH0yNEkX49zJpeV8Hg0rQJSxdMrRKoTzuVCxBIOHCFKC54w8vYIa+LhyqQJm7w+cuDIbUTwEc71uX0YrVyoCicXktYuFLuwfVlyERHK5F6S00j+JjluGILJtzUJxVQqG549cvgQqRQx8RgFMeSQluU+CZHuFOezMtwmnMuFiGWNUSNEORJ+KGUaVHcj5PD0IVbwTrHiC0YVL7l+FOzgDcQrSG4fQphCnsLuDf+jFyAeQnIF8dpRiY4HLN2EXWGGZ4+KeAqByH/ky0FFhL6rAnzJYgyRS8qqMIxTbxRkca/LklLBGKUwifQG8iyjKwmJksP95POEFl10cCAVhNxfvvyQSsIXvvAzI5zLhYhljVEjhCgR5Mrx4UxBBJMeH8hU8hKyJkeR1UMIsVJ1SoiXNjT8j8mRAgQ8hIR5fTEHIodWQ1SX4/lC4JHLli5h2ZIYXj2uDQHL73hHEcKVWdiSKSDwGGvkbPpm1xiFMIxHiosQ1nQYKA3sn9xbWjPRJgkPtigbFFTR1ol7iWjmM4XVnIgUhF0WIJzLhYglHDhCiLKDt468PgQf7XAI85KDhvCjmTd9HCnmoLiDvD68fXhYaF/DCiSIOvLXqExFEGWSB89btJkyQo8cQ3pg0rIn7LspSgZjiZ6JNFgnPE5LGwyBR8ifXpb09CTlwPcnLQm0I6IROv09Wb+b8S3KBkLOf+FhtR68euT90vA/uo67J5zLhYglHDhCiHWDD28xIdKfkG/trJSCxwWxh5eFvnI0xyaU6T0woTjKNPNiD68lwpZ7gehdV4UsVQHGGUvZ0SAbbxAhVQwBzUo4VIrTqJptGIslhQpp1i7G40zol2KlogjDjesSjlWc4xVnm7JQ3PMA8plZYAARTnoHnn7SNVhQIO41CudyIWIJB44QouIgdwpBg8ih/QwrvZCn5sNroSDKZEPs4c1jkiPMjVezoFVdROEMHz7ceYDwkhZ0DykWQqjRwzAKoULaCNHonibuLFnpm8iTK8bfiHHWvWYfrBzkBTnb+ebpPJ/9+NVzgHQEVgdiZR9yBUldoC0RuYOAIML7zWOsZkTLIP6OLplJnhrXR1oDYpZtORZfnuK8kjxGIRTXglji/Pz1eAGGSMXjTg4thVXkw1L5XRyBxnuYa6GROufBOZN2EOaX0nCfMc3KQJxLz549XfsmD830ffU5njxy8/D6s/IR1xgSzuVCxBIOHCFExYEXhFVcyNkjVEPlbCh+qpoRosarRwELXqeCVnYRBYMIwovHmIrLuyMUSHEQ4Vy8RR7ECn0YWQKRLxy77bab8y6xEo7vKUkIkTWS11tvPedl4m/yyfBUI8zwRBMK5rmcAyvLED4GmqS3aNHCCUv6T9I8nbWOWXmIED0CjPcA4eRdd93VeRsRqj6Ej3DCI8lzMb4MsDwk+2CJRlZXQQR6ccZ5IUTr1q2bf05cDxXreMoAwYf4Y9UYtuG4nDfpFAjMuObxHs6HFZjwwnOunAP3jZArq7z4tANELKkIhLc5B45x3HHHuVWC/HrohGsl8kS5Ew4cIUTFwAc/bVf8mq/89CtYhMKnqpi/doRe1KNXlkrQqgieJNYFZ21pvGmIFTxuPI7HiNxQBA1FEwhpQFSwDjFiBUGCeME7xuuAV5BlHRFEhBIR4XjyCK3TIxLBhAeM8C9jmX3ioWbZQHL3EJJUiCIU+X2zzTbLf5zQMe8FjofwYZ1sRBitdUhpQKAhfvB8Id7YL+dYvXp1d3x6SNKj0ws58tcohsIoWkB0MY5YgpJKeLzFPJ/lNAFvIqKWx/HCcVzEIs8j3xDPYxyIPxqjc24s50jqBc+lkTrrdSMm8TYibBHOCFeWfGTZRM6Ze8094EsesMykr6KXyBPlRjhwhBDrHi/wfHjS974LRU9VNR+6JXRNVTKTZ0FL5Ym1QXgdffTR9p///MeOOuoo56VCPCGiTj75ZCdAuL8IIiq8EWD0d8TDhMiIimpEIgKKfFFEC+AdwyNFGBRYy5xKW8RNtAiDfSOC8LYhbvD23XXXXbb55pu715lVXPzx8Qgi3Lzo8eDtwuOHYATEIH8feeSRThh5EKZ49BBUNBBnSULuAdtxfh7EHn0n+cn7kBY+iEP+9iBmKVRB/BVUOYxgRYTxXMSqD4tz7xB6eAIRiAg6PHjch+ga6YS9uV8IOx96JhdXIk+UK+HAEUKsW5jYqKglDIVHhNCOBN7aRk6iX92DiZsJPcxzEvEg8hA3NWrUcPePe4moq127NuLAiQ6Es4eWKHi78M7R1oZiDUQfIpAvI4QvCe327t3bbU/IFU9Xly5d3N/ko22//fbOc0grHEQJz2UfiErOA+8Y+XGEievUqWMtW7bMPz75b4gbn+eHZw5PG541KoA33nhj9x4hdM95c2zawEQhDxBRR7/KmTNn5p83+4t+QSCHDrGFKCP37vDDD3eePYpJOB7nzr6oRMbjeOONN65xHA+ha8Qy9xnBFgUB6HsRUi1fs2ZNJ+DIxWP/HIfzw0NKyJvzIFTOeJfIE+VKOHCEEOsGQk188NNfjAnLC7xQ3Mj+Me4TooMJkubS4TJoIh5y8vCenXLKKc5bRe9BBBLi6ZhjjnHCgTClzxlDfLFO81ZbbeW8eXjlaMGCiDn++OOdJw8BghBhHCOcot4vBN2WW27piofq16/vjstzvQhiH4RzyQ9EhOFVJPzrwdOGCMXziOji+Hgd/TkQ1kTs4WHDk0foE49gFAoaEHmEcBF5zZs3d54y8gQL8gJzDgi8bbfd1j2X43KuGOdIyJsvZHHgLcRLSg4iYfCCQDxzXYSyuR5ELz/ZP9eKhw+vI6+ZRJ4od8KBI4RYNzDB4SFAvCDwMm0Fi3Vh3B+/FBp5elQvpvO6yBWFL7xgnOGl85CnhneNsCahSPIdAVFEjhghWMKmrIxB7h3igxw+cs/wDvJ8xOKjjz7qtvUij6pUvFW8Xn5JP55LxSoeWEKRiBa8XxTUINyiwoV8QcKWFC7g5UJM8nzainAeGEIULxznHOfJozgEoYankOvBc4jQJVcu2o6Ha6DilwIN3o8IWJbZI8zMcfHKUdzBNSM+2VcciDwEM2ItGjYGxChVyoSn8QQSXkZ0ci+5Vn9feR7HZExzPOXkiXInHDhCiPKHiZHcJqpoycHzKxGEokYWbwg9vCxMmCTea1WMwvGePHIaSQ+IggDBC0Z1LGFcQpe0ECEEu9122+XnvnkQQwhDiiEQPNx7tomKPAot8FRxvMmTJ6/xfHLS6KeHSETk0QgckecLHwBBQw4cYgjxE4Xzo5ekLxAhRIwXsSCRx7UhCPEuck4UNEVzDMl/Q9ASLm7Tpo3zIOJtDMUcQozcOgo+4uA+Mi7xKuJFjEL4lSIqvI6EZfFE3nHHHWtsw2cCVb2cN8dG2OGx9s2QuScIbsS4RJ4oNeHAEUKUP3gGmJSYBCXwSmd48vB0kONE2Fb5eQWDtwghhJcprjqUliaEGQnPEnpFsCAsCLVSKYoAQqAgqPG88RgCGw8YVaV4y8izQ8wAQpLKUUQiX2T4QkP+HuKHcCjhT/aPyEP8IFwQN9H2JOSrsS2iC88gbVPw6vFewcPnc/gQjeQHhrlyVOIiuBBwiDzEJiKWEDJrSyMOOS+KLng+uYOsSoOnDy8kBRg//vijOy7FFoRSuX+s7hEHIo1cRDyBGG1gEISIV+4tnlS8goRyKaziPKgk5jxor4KXbscdd3THxeuHt5LjUSRDOxUqlineqFWrlrs3IeFcLkQs4cARQpQf5C8hSAhvkXuDwFMeXumMyZ57yO+E2aKNZMWa4DVDLCByop413z8OgUI+Gon/3FMvBPHYUTyBqEOYkA/pG3QTeuX55PE988wzbjtfeAGEaBnfFHdwbF4vxA/HIA8VCEtSicr/8bRFhToCkvcJXi8EFmKQ3DWEG3l2CDKgSpfWL2wbhTGBhxCPme89RwiWfSAeCflTIcz5IFzpoQcIXoosyEVEjOL5Y/8cl3y8wlbzIPRLrijbcs6MTX5yb/B8+ucihsnfYzufesB1cjyuBwjdcv+4bkQ64WDfWw8xGRLO5ULEEg4cIUT5wQc1LSAQdnihQuEiK5khPPD0MEnj+RDxkAuGsCA/riCR4qtL8RKxAgUQ1sSLRHiUAgnEFXmQCCEfIucnggQvW3Q1C7xyePQIQfI8no8YJOzpz8E3PMajhrCMrijB74QoCaESDqZgASFHbzvCmeyfbTgmBSRhWJcwLF43jucLLRCziE+8dYhLvI142shT9MdGaBISxnvI/zkunj/f8sWfW0FwH9n2sccec8fgJ+I3mjuKGKRXIWFr7gtefbyFnJtv9o3njtcLLx+fG/61oNI5ruAonMuFiCUcOEKI8oGJhg9sqhbJw5MHr+zm++eRu0T4Ltp3TPwDgojxh1AoTKAgMBAyCJJo6BSxhbCjMjcsdGF/PI9ihrjcSI6LgMObGC36AI7Bczkvnht3boguqtARQIQxKbaIwvM4dnQlFPbDduyXx8P9co1cD2KpoMbaCDG8exw3bnm0omAscs38LGiVDH8eCFL6A3o4X67LN3HmHrAPRDHny7WF1xTO5ULEssaoEUKUG3yQ03EfQYIwUR5e2Y17iDePe8rSWFRahpOfSB/K67Urr/1UBpx7cc4/3Cacy4WIZY1RI4QoF/B+UBFHThMCLxQrsrIZbSbwjhJWJJQXToBCZDrhXC5ELOHAEUKUHarqSACnT5mWLCt/I/RNwQBri5Lf5Zv6ClFVCOdyIWIJB44QomyQa0N/LQotqKJTLt66Me4tlYrPPvtsfoK8EFWFcC4XIpZw4AghSg/J01T33Xrrra5dQyhMZOVr3GPuNZWLcUUAQmQq4VwuRCzhwBFClB6qEulw71t9hKJEVn7m++bRfJeeb6zMIERVIZzLhYglHDhFQYIzoRG6l1Nurq7zQuTCe4HebXTiV8uUijGaxyKoqWL2zXJF+eKrP31bD5EahHO5ELGEA6co6NfDUjd0LKcbt286KURVh55g77zzTn7jYxVbrHvjXuMxvfbaa92KDGLdwLJnw4YNc55qkRqEc7kQsYQDpyhoNElriC233NIlPhfUTR2K6v9T1P+jFHe7OIrz3Oi5FLV9Uf/3FHc7kRkwCdLNnrYpWp+24oyeedxvln4KG/eK8mH48OFuFQrGuEgNwrlciFjCgVMUiDyWWdl5553den9FiTwPS7awUDMNYqMfxGxD+4OoRxCPyKhRo1y384KEF8elGzv78x3M2W9cKwUfauDcfZd01r3kgyu6rmMIXktC0nyw0aU82qE8Co9zHuwv2uE9es6co98XS/8oSTyzYNy1a9fO9W9TLl7FGUKaClvspZdeKvT9LEoHn5ksf8Y6uIMHDw7/LSqJcC4XIpZw4BQFb3jWEiyOyAOWnyGMwhqGTZs2dZ6O999/3wkdn89HM1PWK6S3WK9evVz7CbZl/cMOHTq4tQE9iCME1bvvvuvW/7vrrruc6GT5KNYeLOhDHgFIHiHnzvqArNfIGoUPPvigOwbhCI9fg5EFtPk/23E+nBdC1Qs0LwI5F86VbVjzsGfPnvlJ4OyLxbJZ5/H+++932xDmxuvAGowFLe0j0gvywSi4IEyLJy9dvHgXebuorHZRpRjXwP1G5PHZopBt+cNnHOsvsyZr3BqqonII53IhYgkHTlGUVOR9+umndvDBB9v+++/vQiq0PNh9991dqNe7/hFfTI6HHXaYHXfcce53quZ23XVX22uvvdxCzp5ffvnFVdPttttudsYZZ7hQzYknnui23W+//ZzYiwMPG9/09913X9tjjz3cuVx99dV2wgkn2D777ON6bfk1EvEgXn/99e4azzzzTJfvc/LJJ1vNmjVdQj0eRPj111/ddgcddJCbZLgf/H7ggQfaF198kX9cBG7dunXtlFNOcfs6/vjj3T0gWRyBW1okDlMHxt1NN920Rm88vHprGY+vA2sc/B79O7To/93vF16Y2ta4caGGyONzAE8Ty5yJ8qVLly7uyy6f+9Evw6JyCedyIWIJB05RlETkMfEhohA0r7zyivPSdezY0a0EsO222zrv16xZs5xYQuCtv/76bp90sO/atasTZdtvv737APcLUuO9Q+AhrtgX2+EZ22KLLWzjjTd2HsE4EFvNmze39dZbzx0LLx3f+vGo1alTxw455BDn0eN6/ve//9kRRxzhJmw+4OjBxbkzoZCLiMcQDx3ic6eddnJd93v37u2u4+WXX3aCkHPnfClSQfix8sGPP/7oPIHffPONE4X09+I5pQWRR/iE8+W6Bw4c6DxKeBHDRb3FugVxcdVVV7kvMXwRYOysYeeck/vzvIbW4PwL3E//e4OGpbTzL7DzGzayhhc0soaNGlvDRvwsmV3QqLE1atR4bWFVXGuc+xNvnvMO8rME5j2Bl1xysV126SVr22WX2WVNLs+zJrl/B8Zx8ZzyOcAXRlG+EPW4+eabXQRF3RRSh3AuFyKWcOAURUlEHqHUHXfc0SXsRj8ccPnjzUJY9ejRw+XfUa17wAEHWPv27fO341hMmvXr13ceL8Kbhx56qB1zzDEunOpBKF555ZVWq1atAkUeIdPHH3/ctttuOyc4oyDIEKJ4YtgOcYpw4kONXCv2n52d7URp9erVnWeOMOtbb71lW2+9tfMgtG7d2l0X1WcIObyUCC0mf7yRJ510kn3wwQfuf/Pnz3f7RqDNmDFjjXMpCYg8RCTXhVfwqaeectfy5ptv2ocffuiugdeKc2M7jocXkutRVXT5wmtPQRJrqb7wwgvudVjDXn7FXnrlVXvp+efsxcfvtxceaWYvPHqPvfLEffb6Uw+Uzlo8YG8+9aC98+xD9u5zD9u7zz1k70SMvwuyf7Z52D548VH75OXHnX1cYnvMPn31Cfvi9Rb25Zst7Ms3Smb/y7Ov337avn332X/snWftWx577wX79rP37NvvvrdvsG++Wcv4wsaXNcZ3QXmzonTw+ccXV1JW1IcwtQjnciFiCQdOUURFHsKqsGo2wlaEQvGGRUXe77//bs2aNbO9997bPvnkEyeIEHlsH82pYd94RRBITKKDBg1yYV+OG0LOCIKSc4sDYfPoo4+68yGHz+OF0lFHHeWO5ZdHQlQi9hCEeADJ90GE4oFkEsezyHk3adLECUfCwHj6Xn31Vfc4YQ32jTAll2XPPfd094wQMx4H9k0xCvelLGFXrmvAgAH2ww8/OM/iE0884e7tbbfd5vq1ER7G+BuxzQc238wRwzwPgY3gJO+GMDX7oyiFey9vYMlgQkTAI/T5uYYt+tPmz51jc7OH2JxBP9ic/l/bnJ++tfk/f29//vJj0TY05rFffrS/hraypcPa2LJhbW3Z8Dzj9+jfoUX/P6ytrRzZznJGdXC2usTW3hKjO5qN62SW3dlsXCltYhezSd3+sYndzCbkPTbnV7MVEm8VDZ9L48ePtxYtWrh0Fj7zROoQzuVCxBIOnKIIRV5Bb3zCmYSnyJOjoCIqZAidInQQRnicCIeedtppLuxCM1kP3xwRRYQ/EXl9+/Z1Ig3REvLee++50GlUwEVBtCG28B5GvX2cF8tQHXvssS4PkGNOnz7dFUggzAjt8jheSc6DYzz33HNu4kYEcV4k3HP+tWvXdt5EcvJatmzpJntEHD/JTUQEImzx/pFDeO+99+bn95UGzp37jLjgdeDc8UQi1PAQ4n1EcHJ9ePO4R5z7fffd58Iv5CTijeV+EkLHG4jXCeFN6BrPCOFfvCO87hyHa8aLWVZxWiVZPNtsSl+zbARRp1wRM76YVti2/n/8jG7n/w4tfN74zmtvUxILz6c8Lbtj7s/fh5qtKjhqIMof3ut8eSRSQIRFX/pSi3AuFyKWcOAURbSFCnlxheVoEMakyIHK0mjLELxbCAwEIOEWL/JIVsdb58GrRriWwgqqZmlRgoCiaMK3QgE+fAhV4snjQykOhCUiZpdddllrG3LlEHOcE+f2wAMPuBArf9M6gGpZxBjfaCm+oILWCyq+6fI7BSFsyzYIUfZHziBhUYQS14LHDPFEmJfr5Vz4hryuQQhSoYzw4zo4VwQz+YycI9dPqJFw8uuvv+5CjngE8WDec889zjOIMMQTiEik4S9Cn32wPzySCv8WwR/TzSb3NRvT3mxcx1yvV77nq1PZbGzZLDm2Y9lsTEdLlNVGd4g3PIUjW1tiZFtLTP3ZbMmC8M6KdcTy5cvd5xnGZ4e+1KUW4VwuRCzhwCkKxBX5L3i08EyRG4fYQUQRhsVjRrgPYYGYQSxRODF69Oj8fZArRg4clav0n8NjRFUt/cUIIXoQRlSkHnnkkS6UyN8ck8ILPGUcF/GCUKR4Ag8ZvcriwJvGN1IKJwhjcq6AQCHkioeNMCZiEoFGWDaaH4h3j6rezTbbzOXk4dVCvHJt0e24D5zzNtts4/5PIjjJ4Qij6DdhBNRGG23kziUV4HXlXiDaENqIP14nQsB49xB3iFvuIR5Rfj755JPuXuDJRBziKcQLiGCkeIZiFAQ89xQhXGWXREosNJvzs9mUbmbTe5nN7JNrM/J+ygq3Gb3NpvXItT/G0EgpvMNiHcDnLcVhpKdU2fduChPO5ULEEg6cooiGa/FE4elh0sdLhgeM4gTEHUIPAUjlKduRE8bEjwgg941QLcIHEHZHH320qzjt379//rH4kCFEiugiZwzwwrEtFascl951eACprqUoAnESByIPQVKjRg0XhsVbhThDrPA31Xl48RAjhDIJvSLAyJ3Dk4dni3xAqn1pJ0DOHQnfXBv5fFwLnjE8c4SEOSfy3RCwFIpwveyDbajOxZPH8xBS6QIhYcQ2nkAKShDFiD9ee0LP5C1y73yPQx7nfjAeuE68mHj/uCeEkRH+hLvxEjBe8Aams7fA50eSntC3Xz/r17+/9e3a3np//Yb1+eRp6/Pps9bn8xesz+fPy0pj3MPPnrM+rT6xfr17Wr8BP1m/fv3d/cbCRuui9BChIXJCwQUpNSL1COdyIWIJB05REHalnQhhVAQPYgujFx6FCXjoCKciypiwmcwRePwPQYQXjkIKPER44gBvD3lhCEY+WDx46RAN5I35JscIAcQZRRpU2lJ5i8DAw0cBBEIyDrxUeKLwQHJ8vIacD55Gws60H/E5boiY6667zl0L4gwxxjda8givueYaJ2bJd0P08AGIF5KwNNdGqJYwNQUk3CtEMcIOAcv+MI6JeEUkpmOY018XEyr5eghjRBr3GGHPa0jVNNdHpS/ilpxG7jNtRriHCGgeQ6jjBWRbBDXVx3iDEdHsn3vMPeJ4eEJT1aPAeRHuJmVgq2rVbJuaNW1bbJsatm21LaxmdWwrWZlsS6tZbQvbtvpWztO+ba1aVrNmLdtqq2rOc857VstulQ98wcYj74vEROoRzuVCxBIOnKJACPEBgBBi5Qg8OnjP/E/CpST6+2/UbM/EjxcHzw8eHQRVtBwfMYfXi/y26Ddxv7oF/2OiZzvCnyNHjnQiy3vZ8AghvBB5BYVrCaPitUOIUWGKxwWvIvloHDeaW8jvPEaxAtvguWTy8MUMHN83BeV8uR5ECoKPbfl/dH9cB8KH/3EP+PBEBEbvUabBMm4Up+AdJQSMkOc+4cmj2pr8RTyaiEC8gVQ+I9YpAOFLAQKQFUJ4zXi9eB14Hq8/4ym6CkoqwOtN7hJvKVnlGLm7jC9RdsgjJmLBl3Gf2iJSi3AuFyKWcOCkMniNCAniKUTcea8OwtGvllFQM1Q+qPAUEtIlV64gylNwsa/i7K8422QavHYIZTy0NITm9UT4UYXMMnGIP3IBCbEjnvAGeqPIhjwhPIDkZrJMHnmZfOlAPPPFAFFP4nhFwfV8/e13dvQxx9jhB+1nRx64lx150N65dnDeT1nZjXuZfz/3cr8ffsiBduRRx9gDDzZ3XwhF2eFLGd0TiIwUVlwnKo9wLhcilnDgpDKE6wjvEe4kp43JHi8Q4g6Rh7enoGaohBLxCvHcwkSeSB0QaYg1inPw0CDkEHXk+RFGouoXzx9G9S9jA1GI+MNjSoEQoX08vniXac+DRxZvrA8Jl1dbCGT65KlTrd33n9sPbz9prd96xFq/85i1fudx2TqzR909/uH9F+zH77+2IUN/sz9TzMObjvCeIMWE7gLRQjiRWoRzuRCxhAMn1WFipoUL3zIJz2CsVEE4ONpWJYSVOdgGoUeoWaQPeDrxJlDRjPAjdE+uHuFuUgfwBlLMQcoAIWDCu7R7wevr8wCZsMjvJCRMbiaVwITqEX+E1wkB80WAEDP75MsCx+GYxcoDxGu7bIHlzBxiORO6Ws74zrk/ZevWxne0nMm9LGdhtiVWLsn1noevjSgRpKXw/uCLlDyjqUs4lwsRSzhw0gEmerwxeGYwJnkmY08Y/vQigTwu39Mt3EakN3gfEPJUUdPuhqpdcjXJhSSEjyeQJtiIQFZHwQtMzh8FNRSBUBGMKKQqmMmNlU7InaSABBFIkVC012MYik/+Pc9sxuDclRomdTWb1F1WIcbqGF3NpvYxW5BttjK+ObsoPnQ44D3w+eef5xfHidQjnMuFiCUcOKlMcYRZOPn6x0LiHhOZDZ5e38B6yJAhbjk7WvKQ00monzY4tOQh3w8R6HsCUhRCLiC/03qHAhK+aPgRlFwy3xJTBpmN62A2uZvZ1F5m03rl/owaj8U9Xhwr7fNKYnHHiHusOFaSay3udqHlH6OH2eSuZpN7ms0dZ8mV6VexnkoQKaH1Fbmy5ZXOIMqfcC4XIpZw4AhRVUH4EwKmSIfcPap5EYB4/ij+wMtH8U637t1t6fK81IBliyw57We3IsOy336wBb/8aHMH/2Dzhvxg80P7Jc/Cx4tjpX1eSSzuGHGPFcdKcq3F3S60vGPMHfy9zR38nf015BtbNaazJedlm60uOHVDFAzeakK1iDzyYUXqEs7lQsQSDpySQqsMLFqBVRW8ZKHHMPxbpD/k4jGumfjwaOC9QwQuW7HCcnitly6y5PQhZmM7Ws7INjbku9fs5fuvs2ebXW3P33utvSCrEHv2rqvs2buvtlZvPGyz+nxqNqW32YJJEnqlgNxUvNV4raNtrkTqEc7lQsQSDpzigqChcS15GySw0wuN3krkQzE5IvzIi/J/ZyJ8CHKNeH4KK/oQmUdy2R+WmDbEEqxFO76jrR7dwd586GarseV/rdoW/7XqsgqzaltsZltvtbk1OedEG/zjG2bTuptN6GGJ+RMl9EoAn+m0TGGlGkK2fIaL1CWcy4WIJRw4xQHPBm1IaGOy33772cUXX2yvvvqqW9uVqkYED20vqGak3YUXefReIgRG0UQm9F7iWqjspaebGoZWIZb/aYkZv1rOyLaWHNXWJf6vHtPRHmt66VoNemUVZ8fW29/6fvOq2czeZqPbWs64bpZcMMmSOf8UZYmC4XOadaj9ij6Z+uU8UwjnciFiCQdOcaB9RYMGDaxWrVpuiSq/LBUVijSwJaxFc1uW+mI5Lw8J7CxRlim9l2jCu+uuu7r2HHwoiirA8r8sMfM3yxnd3nJGtLbEmI6WHN/FcsZ0tPbvPG6XnH28NTrtGGt8ev18a3TaP7+X9W9+D/9f3G3i/h+3bfiY/7uw5xT1/OJsU9TfBf/vGLvojPr28v3X24Tun5hN7mHJ0R0sZ3gry8nuYcmFk81yVEBQFHju6DfJ53RFNhIXpSOcy4WIJRw4BeHzzXjz05OMtWhPOeUUF7L1PcvI56CFBflLtJ6oU6eO8+7xjZD2JazXiueP5cLodxf15vF8QruEfOlTFj2u75Hm1y/l/3jOOC77LigXjl5nLKvFOrpsH/3g4nk+zwr4gJs6darbPlzk3O+f57AIPddJuw5WZ6hbt65bxJuwrchkkmYr/rbkzGGWg4AY0Sr35+gOlhjTwRKjO9iSX1vZvJ++tjn9v3I/vc0d8M/vZf2b38P/F3ebuP/HbRs+5v8u7DlFPb842xT1d8H/+8rZH4O/sxUj2llyTMfc12ZUeyf0EuN7WHLRVAm9QuCzjc8wlhPkS7tIfcK5XIhYwoFTFOTe0VCW5cFYPYJ2E4g+BBItKWgvgRjr1auXE0A0o0WM8T8E3mabbeYqt1iCiseBthYIJjyBiEK8gfQ48+ECRBfVjjQzphM7x7zzzjudkIzLG0GUcT6su4gAa9KkiduetWO9gCPkzDFYSouGn+ybfmmcL2vQsiJCFM6FxrlUV7INnktWVzjyyCPdigsISZHBrPjbErNG5IqHka1zBUSeyPNm4zqbTeiaZ91kFWZ59zy7s/OsrvG6jGznPK45E3pa8s8ZZon0TxNZF/AlnDy85s2bu89CkfqEc7kQsYQDpyjonXTppZc6sbbNNts4YcaHQseOHe3CCy90AhDhRUNNRB7CiQ8QlqTac889bb311rNzzjnHiSS8czSaZZsDDzzQzj77bLfvI444ws4880wX8gU8dwi2+vXrO1HFtkcddZQLl8ZVgCEoOReOz09Cyoceeqh7Hn3QEJWINs6B/Z166qlOXDZs2NAtkbb//vu75GM8ix7aChx++OHu3C677DJ3nnvvvbftvPPOTvjhARSZCX3XkrNGuZCsEwyj2q0l8HI9eh3NxnaSVZLle/BCG9nGCfPEhN6W/HOmWVK5ZiGscoHAY7UYviCL1Cecy4WIJRw4BeFDlvRO+vjjj2333Xd3Qqtdu3YuH41CDERco0aN3IoSLBSPyEJg4dnDa3bsscfabrvt5jx1eM8I2dJZndw9hBJeOlYnwEPmxRTbsT9EFwIRocaHEYvYDxo0KN8b6OFYJA5zLghOvH2ElAlB1KtXz1mfPn3ctjy2wQYbuLw6mt9yLR988IELRR900EGuatgLUQTejjvu6NZOZS1URB+eyY022sj1Tws9fyJDWLXMknPGWM6YzrkCj1y8UETIUt9Gts0NsU/sa8m/ZuWG30U+FMqxVCR9IaMru4jUJZzLhYglHDgFEc17Iyft6KOPdt4vvHZ4xQh7Hnzwwa4Igce8yGMbD1W4xx9/vFttgP2xhiziCc8Z2xPyJC+EwozrrrvOttxyS7cKAR86Tz75pG277bZOABZW9YXoQ4Ti+WMdUvLy8Mj17t3bFX3ssMMOLmwLNLnFG4nAXLBggXuMHD2Wu0LAERYmx5AKYUQjgjQK57L11lu7pbHkyctAVi6z5NxxlhiLwPsnB0+WpoYHdkQrS0zqZ8nFc+TRy4PPV774ssYzaTQiPQjnciFiCQdOccArh6cN0UYolcII1gU95JBDChV5F110kQuP8j+KHsgB2WeffZxQ4gOGlivk6+H9Y/+bb765a0+CUGNReXIAv/jiizXOJYoXomyPYKTFCWIPzx8eRkKr22+/fX7OCZ68nXbayQlJD57A5557znn8aCfAuriXXHKJHXPMMa5gJIoP9yLyVHiRYeDBm5ttOQi8AnLwZGlo5OjxeiL0/kboyaNHMRlfjPlyyxdwkR6Ec7kQsYQDpzjQ7+6www6zxo0b25w5c9w3weKIPLbneRQwsM1HH33kRB4hXMTdbbfd5jx45Ogh+Jo1a+aKOvCy4cnD4/ftt9+ucS5xsBg9xRaEdsm3o/iCxefJt9tll11cw09A5NWuXdt9wHkozMBDh4eR3n+sc3rWWWfZCSec4HL9orDuKRXDiDy2ExnC6hWWnDfecsZ2KTQHT5amNoIcvba5Qm/xvPDVr3LwGUs6DKkq6veZPoRzuRCxhAOnOFDp6kUeRQzFFXl40/CI8W0RTx6eMIoXTjrpJBe6RTwSqqW44+eff3b7JNSKt5DF4gmhUsBREHjyCCVTFIF4QyRSnME+McQYghJxBl7kUU3rQeSRn8f1vfLKKy4JmWsgFM3qHlHIQ6Sgg/2qujZDcAJvgmuk6/qsSeBlpvlijCkDLblkfjgKqhRENvhyzWckOcgiPQjnciFiCQdOcUCM0QPv/PPPz/fkEcqkYIG8N0QeHxgIIwSX54ILLnDetfbt27sPk+HDh7tiDDx0FFx4yKX75JNPXEEDgpDqXFbSoEDjyy+/zN8uhFw9cu/InyNnMPqBRa4JxyI07MOueOoI3yIEPYg8vIZcHx49qnf5hss+8Q763D2EJ7l87A9PIeJSpDk5Ky0xf5Il8gWeQrSZa+0tZxTFGHlCb+k/vTmrEnxm8pl67bXXutQUkT6Ec7kQsYQDpziEIg+vnBd5eL2ohqWFCiKPfDYPHq/tttvOfWvEc0f+G14zCirwChJGpQ8f+XdU75544omu7x6NhxF5PFaYyAOqcQnR4vWjOpYwKknFhIC32GIL++9//1tskYf30DdAZuUOjk+LAbyUhJrZplq1aq4Hnwov0pycVa5hbk5299wiC3nwqoY5j14bS0wdZLbsz3BUZDx8tpH3TFEZv4v0IZzLhYglHDjFgXJ7BB4eLLxuePJoK0L/O5bFoUIVzxzCCG+cr4bFm0ZuW82aNZ3go+0ILVmoXqXQAS8f4U/alxDa7dSpkyuiwGuGIEO84QUsDLyItAEgV479kIfHedBOhfAthRIIRc6JIg7+H83zQ6C+9dZb7vqind/J86M4hBxCzpH/n3vuuc44f/pMiTQlmbDEH9MtZ3zP3ElfAq8KWe7ydBRkJKYOtuTyqiX0aDHFl1o+E4mYiPQhnMuFiCUcOMUB4UWjYvrUIei8twuhR3gV0UdYkw8QknrJlcPwktELD28Ywor2JMC2/E1VK94zKmKj4VuEFxW9VMsWR0zRfw9BSQsU9omnjjw/8vsQjohUzpnwBH9HQ62c+9ixY931ccxo6xhyERGbeB/5UMRDiMeSPnxcW0HLq4kUhmXx/phhORN7q8iiylpE6M0YalaFhB6ftfQV5bNO/fHSi3AuFyKWcOCUB3FiJ+4xKOhxT1H/X9eU5Pgl2VakAMmkJf6YZQkn8BSirdJG/qXz4ra1xMxfzVasuX51JsIXXb5UI/LCpvIi9QnnciFiCQeOEFWCZMKSi+e6FRDkwZPlWp7QG93ekjOHma1aGo6ajIEvpLRLoTce6SYi/QjnciFiCQeOEFUBGuHSJ831S8PWmvBlVdLw6Lnl6zpaYuZvZisz08NF7jLrjT/44IPF6j0qUo9wLhcilnDgCJHp4MFLTO6fu/KB89zETPayqmv5Qq+DJX8fkZGhW4rZaBFFJwNylEX6Ec7lQsQSDhwhMprF8y1n8oC8EK08eLICLF/otbfErJEZ59GjgI22UizpiFdPpB/hXC5ELOHAESIzSVpyyUJLTv4pr5JSAk9WhEU9erNGudVQMgEKLuh1yupE33//ffhvkSaEc7kQsYQDp7woTqWpb61S1GNClJXkskWWM2Vg3ioHCtHKimlO6LWynLGdLTFnjNmq9Bd6tKxieci77rrLevXqFf5bpAnhXC5ELOHAKQv0z6Nii585OTkuDED/PFbF4G8Pfezooffdd9+5FS6Abfv16+eSgOlpF91eiLKQXLrIEtN+zhV45OFpuTJZSYzKayf0ulhybrYlVy8Ph1haQU88WqfQ7H3SpEnhv0WaEM7lQsQSDpyywHJhzZs3dw2G6btEI2TWd6V5MGvAQnZ2tj388MNWv359t2wZK1ggDGlczGoYrDnLihMslSZEWUku+9MS037JX77KtckIJ3GZrCjLW+c2Z1wXS8wbb5bGQo8m7iwtScN45eOlL+FcLkQs4cApC6wvu9tuuzlRR/UWS4Fddtlldvvtt7sQAbz77rtWvXp1O+SQQ+yGG25wXrsOHTrYDjvsYPvuu69dccUV9sMPP4S7FqLkrPjbEtMReHjwJPBkZTTn0UPodbXkvAlmOSvDEZfykArTpk0ba9KkiQ0fPjz8t0gjwrlciFjCgVMSSOBlyTG8diyJg7hDvNFcE5HHt0Q8dywbxv9Xrlxpd955p22//fb22WefubVrWS/xlVdesW233dYtd0bV199//+32TciWn9EcPf94mLvH8mosLca5FOYFXLFiRf520ZAw+/L7Bh9uZnt/HH+9PJ//KXcwhVmx2C1RFa2SXGvSlslKaj5HL7u7JeZPNMsp+LMmFeHL9nvvvWe33XabTZw4Mfy3SCPCuVyIWMKBU1wQbKxDi2jDU/fjjz+67unHHXecE2uEZ7GffvrJfWMkD69Vq1Z28MEHW7Vq1eyOO+5wyb8fffSRnX322bb55pvbpZde6kK948ePd2vXsv+pU6fmHxOh2Lt3b+f5I8TrGTlypOv3RAiiadOm9umnnzoBGUL+H+dGVdmtt97qzjkariB/kGNyvvyPUPObb75p8+fPd//nf3fffbd7PmFpPJXKHUw9ksv/suSM33InZdcqRQJPVo7mQ7fZPSy5cJIl08ijx+fxE0884Vqn+LXDRXoSzuVCxBIOnOKA9w2xhKA7/vjjnTi7/PLL7YADDnAirmXLlq74on///nbNNdfYQw89ZNOnT7ePP/7Ydt55Z9tss83sjDPOsCeffNKFeI888kjbdNNNrV69eu5bZvfu3a1FixYu1Iuw8uCt4xto48aN7eeff3ZCE9F3/fXX2ymnnGIXX3yxXXDBBXbEEUc4sUeCMd43oLs7zyPv7/zzz7fTTz/d/U4CsheEw4YNc2GM8847z0466SQXPuYchgwZ4jyP5557rvsfx+HaOW8ainI/RIpAiPb3Yf+E1kZruTLZOjBSAPKF3pS08eh98sknrqr2l19+cZ+fIn0J53IhYgkHTlEQvsSTtvvuu9s+++zjPjSoikXsbLfddlarVi17++23XciVQoxDDz3UCSOEFJ60s846y4VrCe3yrZLnkptHuBbvHt3XBw0alC8a8fZ5CJUizvbbbz/3PDyFCMwdd9zRiUJCwzyXx04++WQn7AjdIt4QduQLcp78TdLxhRde6ETnBx984AQkYo5w83rrredE42OPPWYffvihKw6pW7euu96ePXs6TyOPc22nnXaaE7AiBVi51BKzRuROwvLgyda1jcz9IpEY38uSi6alvNAj6vDcc8+56AuRFZHehHO5ELGEA6coyK+7+eabba+99nLePL4N8uGBAEOsIbgozce7hcg69thjnZiijQpeNb5FIpa6dOninkcIlv0gGr/44gu3DR44QqKHHXaYC5t6OEbDhg3t8MMPd2KO/BJE1p577umeyzH48KKYg/3THgCRh7eudu3a7vxGjRrlnkeogn3vscceToQiEHneMcccYwcddJATl4RyMdq9cL1UBJPHwvMRmAMGDHBJzD6cKyoRJ/BGWs6YTsrBk1Wc5XmMExN6W3LRdLNE6qZv8GWUzgZ8eZUXL/0J53IhYgkHTlHg5ke4nXrqqS4nLQoFFIRKX331VSeCyK+jTQrhTQQYkOdGGBTh5Hnttdec2Prqq6/c32PHjrVrr73WjjrqKOcN9CDyEIw8TiiYv3nu/vvvb3vvvbedc8459uijj1q3bt2ciPPFGQjGDTbYwHn3yKlr1qyZCyFTyUvoGNFJ2JdrI/x8ySWXOG+fh3PHs7jFFlu4Y7M/qoQ5T0SqDwmrEKNySDqBN9pyxnTO8+ApRCurQHOh21aWM6GPJf+YmZJCjwgMX7offPBB+/zzz5VLnAGEc7kQsYQDpygIsRLSJPeN8GYUknnxrCH28HYh8sh7Q+TNnj3bbUOoAEHG/wBh9OKLL7pQKqFfIGQbJ/Lwql100UXucYQcH1SIyS+//NI9fuCBBzqxSG7fAw884Lx2iDDE3IYbbujyAPkdwUauHWHdK6+80uUF4vUbOHCgE7Cs6chzo4wYMcJ5BBGBhIs5X0LPiD3vyZPIqwQQeHPG/uPBk8CTVYbRZHt4a0tM7GvJP383S6aWiCKyQtcDvgRTWKbPqvQnnMuFiCUcOEVB1SthS4wCiSiIO0KpoScPz1hhIu+ll15y4VqqYgEPGQUb5Ly1bds2f/8INsK1iDg8gYQcCEFwLIQW/fX4pkpBxC677OK8hlTMUphBRS9hChqB4gUktw5r166dW9qHb7o8jshje0Qd4KXDK0hOIcKVSl7ELKFfhB7haVbpUFPRioeVBxJzx7olp5SDJ6t0cx69NpaY1M+Sf80yS+Z6+FMBPiNvueUW13SeFlAi/QnnciFiCQdOUUybNs1VuJLjRnUsIVOEELlwtDBB9CDyqK6lQAMxSFWrF3k8F0FH6MDDBw8FEFTfAnl/tDjZddddnaeMvDpEGCFUvIjkzJGTx7dTzoH2Kewfzx7bkk9Xp04dJzAp9nj++efd+VLYwb45XwQi4Vk8fhSK+GpgvIRXXXVVfqNQ9ocwvPfee524Q8xxHD4oCf1uueWW7hzUjqCCWb3SrTyQM7ZrbqhMOXiyVDCfozd5gCX/nsu32HDkVjh8kaYVFREKWlaJzCCcy4WIJRw4RYHowYtGyJLcOrxvLF+GsMN7ts0227jfEWB44SieaNCggfOoAa1N2I6KVeADiIpXqnKpWAUqc8m1q1mzpiuKoI0KghERufHGGztvIcdEWJIfh7CkRx/eN8KuiC8EIsfiuHjfyNfzwhSxSHUuYo5zQWTiJUTM4T0klPvbb7/lnx/tWgjTEg7++uuvbcKECU4QUoG70047uTV4C2vALMqZ1SssOX+i5YzrrhCtLPVsFEvotbacyQMtuYSVfipX6PF5SiEan5E+giLSn3AuFyKWcOAUB1aLeP/9953IQcTRP46CCDxniDJ6yiHyCIPyOGEC7+lCZNGDjkIHQETRG48w6ffff5//GDlxiDraqBB+bdSokQvhckxCpYg8vHsITo6BOOM8aJWCN45vrQg1X/2LqPSVuWeeeaYrHGG/VArj0QP2Sa4hXjtCxh5CHTRF5vroj8c18nzOhZYE3kupPJd1TzJnlSUXTLZEdp7AI0QWTrIyWWUbIm9kW0tMHWS2bGE4jCsUvvgSzWBNcL6giswgnMuFiCUcOMUFzxXfDsmxIwSLF45viRREINDoO0e+HDlvCDHCukBxA+FU31sOYURolMpa2phEweOGlw1xR5d2hCGGeIz2psODh9gixIv3jkIOvHdREHoUihDaRQByzghVllHz4PWjJQoePb/WrofwLCFmvg0jDAn9kouHAPRI5K1jclZacsEUyxnfMzdEKw+eLGWtY+56yQi96b9Yctk/nxMVDZ97fDHmiy6fyyIzCOdyIWIJB05piBM34WPh3/6xuMeLS/jcgvYXPhb+XZLHoKDHxTokZ7UlF02xxPgeljMcgaccPFkaGA2TR7XLE3qV04CYL9mkoPAFW2QO4VwuRCzhwBEi5UDg/THDcib0Ug6eLP2MlAIn9IZacvnfkYG97r8sko9Ho3iiF76YTGQG4VwuRCzhwBEipUjkCjxWFPB5TmtNojJZqlteix8n9Jb9GY7ydQYpLxS20aZKyy9mFuFcLkQs4cARImVI5LgVBFhJgEazLscpnDxlsnQxJ/TaWmLGULMVf4WjfZ1AMRv5w+Qyq5dnZhHO5ULEEg4cIVKCZMI1lKWxrCuyILcpnDRlsrSy9rne6FFtLfn7CLday7qE3GG6Atx0002umbvILMK5XIhYwoEjROWRl6NEAc3ieZaY1D+vTYoEnixTrF1eykFHS8waaclV6867tnDhQnvqqadcSyhy80RmEc7lQsQSDhwhKpekJZfMs8SUn/KWiWodM1HKZOls7XPTD8Z0ssTs0Wary7+tCS2jWGecdblZqUetUzKPcC4XIpZw4AhRmSSXzLfE1IF5PcYk8GQZam75s1ZO6CXnjDFbVb4ijEbxLNf48MMPu9YpiD6RWYRzuRCxhANHiEqBEO3ShZaY+nOewKPIQr3wZBlseevc5ozpaMm52eXq0aPIgsbtLVq0cGFbkXmEc7kQsYQDR4jKILl0kVsCSm1SZFXKEHoj21hibGdLzh3nVnUpD1hG8vbbb7c33ngj/JfIEMK5XIhYwoEjRIWz/C/LmT7kH89GOBHKZJlsPnQ7rqsl5mWTUBe+Q0oE+XcDBgywZs2auSUkRWYSzuVCxBIOHCEqlOV/ub5hLjTrBJ5CtLIqaAg9lutD6M2f5FZ5KS2zZ892rVOeffZZGzp0aPhvkSGEc7kQsYQDR4iKIrlisSVm/pobniUHb1TM5CeTVRXznuzsrpZcONUsZ1X4likW2dnZdtttt9nnn39uf/5ZcatriIolnMuFiCUcOEJUCCuXWmLmsDyBx9qe8uDJZLk5eq0tkd09V+glSh66HTx4sF1++eXWpUuX8F8igwjnciFiCQeOEOucVXkCz4doJfBksn9sVF5/yPE9LYHQ803CiwFNj3/44Qe77rrrXOsUkbmEc7kQsYQDR4h1Ch68WSNdx3+XbC6BJ5Otbb4R+IRelvxjhlkiEb6TYhk1apTLxaN1CmFbkbmEc7kQsYQDR4h1RXLlMtfh301i8uDJZIUbHj1aCjmhN9MsUXQxRocOHaxp06bWqlUr++uvv8J/iwwinMuFiCUcOEKsE3JWWmL2mNzJyxVZaD1amaxIyytKSkzsbck/EXqF5+h9/PHHdtVVVzmPnshswrlciFjCgSNEuZNYbYm5Y90STq5NhDx4MlkxrX2e0GttORP75gq9AsBz9/LLL9utt95qc+bMCf8tMoxwLhcilnDgCFGu5Kyy5Pzxrv9Xzogf5cGTyUpjeL9HtLbE5AGW/HtubDHGkCFD7LHHHnNCjwIMkdmEc7kQsYQDR4hyY/UKSy6YmCfwyMGTwJPJSm1uTee2ljP5J0sunmuW/KcYI5lM2ocffmgPPvig9e7d21auLJ/l0UTqEs7lQsQSDhwhygXnwfMCj8lJAk8mK7O591JbS0wZuIbQW716tT3yyCP28MMPu3VrReYTzuVCxBIOHCHKTCJhyQWTLSe7e64HTwJPJisn8zl6bS0xdZDZkvnuLbfwjz/t7rvvtueee8559UTmE87lQsQSDpxMRh9+FUBitSUXTLHE+B55ffAk8GSy8jWEHh69Nmazh9nKxfPtl6G/WrN777OPP/oofEcKsQahBhAZTjgACmLBggXWtWtXe//9913uR+fOndeq4Eo1EbVq1SqbOnWqjRgxwmbNmpVy55dxEKJdNNV16ncCD4/DWhOUTCYrFxveymxqb1uQ3dfea/mytXj6Gftp4KDwXSnEGoQaQGQ44QAIIdfjt99+c7kep5xyih111FHOTj31VLvvvvusT58++Um+qSaiaCmAICWM0aZNm/DfojyhTcqi6ZZwAi8vQTyclGQyWfnZqPZmM3rbtL6f2R1XN7L333nLFv6hqlpROKEGEBlOOACiINrGjh1rDRo0sG222cZOPPFEl9zbvHlzO/30023LLbe0I444wnr16uXEYKoxd+5cu/POO61u3br26quvhv8W5UUy4ZZeypnQJ1fg0ak/nJBkMlm5WmJMR7PpvWx0u7fs8nNOsFafvcPX8vDdKcQahBpAZDjhAIhCr6WWLVtazZo1XYPNX3/91X7//XebMWOG8+7RWX2TTTaxSy65xCZNmhQ+fZ1QEm8hIu+OO+6wOnXqSOStKxB4f820nEl9LWekcvBksoqwxOgOZhO72fIxHa3ze0/YtQ1PsR5fv2K2YobZisXhu1SIfEINIDKccABEQSQR6tx8883t+++/D//tRB/i75577llD5E2cONHeeOMNt4bi1Vdf7YTW//73P1u4cKH7/7Jly2z48OHWsWNHl+f32Wef2Y033mg33XSTW3tx6dKlzoP49NNP27XXXmtPPPGE/fLLL+65iDxaBAwYMMDatWvn9oGAu+aaa+z++++3Hj162PLly9228+fPt3vvvdeOPvpoe/PNN/PPj30QZibczP45z2+++cb++OOPNbYRRZDIseSfv+cJPEK0rdeajGQyWflbYnRHsxl9bHLvz+21B2+wx29rYiM7vmP2ez9Lzhputlzr1op4Qg0gMpxwAEQhp+2tt96yrbfe2ho3bmzt27e3MWPG2J9//pm/DQUZEyZMcMIM6Lp+ww03WL169ezQQw+1ww47zHbbbTc78MADXdEGz2W/rLF45pln2sknn2wXXnih25bjnHbaafbMM8/kh4T3339/23nnnZ2Imzx5sjsGxyNsXL9+fbf9eeed58LGbHvGGWdY69atnUhDVCLyyCH0Ii8n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'

