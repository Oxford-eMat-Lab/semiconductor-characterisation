"""
eqe_helper.py
--------------
Helper functions for the External Quantum Efficiency (EQE) teaching notebook.

Scope
-----
These functions implement simplified, 1-D physical models of a crystalline
silicon solar cell (absorption, carrier generation, carrier collection,
reflection, spectral responsivity) that are good enough to reproduce the
*shape* and *trends* seen in real EQE / spectral-responsivity measurements.

They are written for teaching, not for precision solar-cell metrology:
  - the absorption coefficient of silicon is a coarse interpolation of
    literature values (not the full tabulated dataset of Green (2008) /
    Schinke et al. (2015)),
  - the AM1.5G spectrum is a smoothed analytic approximation, not the
    tabulated IEC 60904-3 reference spectrum,
  - carrier transport is 1-D and uses constant (depth-independent)
    material parameters.

References
----------
- Schinke, C. et al., "Analysis of the Quantum Efficiency of Silicon Solar
  Cells" (lab manual), Leibniz Universitaet Hannover / ISFH.
- Bothe, K. et al., "Accuracy of Simplifications for Spectral Responsivity
  Measurements of Solar Cells", IEEE J. Photovolt. 8(2), 2018.
- Quokka3 Modelling Guide (optics: T_ext-Z model, EQE = T_ext * IQE).
"""

import numpy as np
from scipy.integrate import quad

# numpy >= 2.0 renamed trapz -> trapezoid; keep this notebook working on
# either version.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ---------------------------------------------------------------------------
# Physical constants (SI unless noted)
# ---------------------------------------------------------------------------
H = 6.62607015e-34      # Planck constant, J s
C0 = 2.99792458e8       # speed of light, m/s
Q = 1.602176634e-19     # elementary charge, C
KB = 1.380649e-23       # Boltzmann constant, J/K

EG_SI_EV = 1.12          # silicon bandgap at 300 K, eV
LAMBDA_G_NM = H * C0 / (Q * EG_SI_EV) * 1e9  # bandgap wavelength, nm (~1107 nm)


# ---------------------------------------------------------------------------
# Photon energy / wavelength conversions
# ---------------------------------------------------------------------------
def photon_energy_eV(wavelength_nm):
    """Photon energy E = hc/lambda, in eV, for wavelength(s) in nm."""
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    return H * C0 / (Q * wavelength_nm * 1e-9)


def wavelength_from_energy_nm(energy_eV):
    """Inverse of photon_energy_eV: wavelength in nm from energy in eV."""
    energy_eV = np.asarray(energy_eV, dtype=float)
    return H * C0 / (Q * energy_eV) * 1e9


# ---------------------------------------------------------------------------
# Absorption coefficient of crystalline silicon, alpha(lambda)
# ---------------------------------------------------------------------------
# Coarse anchor points (wavelength in nm, alpha in 1/cm), approximate
# literature values for crystalline silicon at 300 K. Good for reproducing
# the correct order of magnitude and the shape of alpha(lambda) across the
# 300-1200 nm range relevant for c-Si solar cells.
_ALPHA_ANCHORS_NM = np.array(
    [250, 300, 350, 400, 450, 500, 550, 600, 650, 700,
     750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200]
)
_ALPHA_ANCHORS_CM1 = np.array(
    [1.5e6, 1.0e6, 5.0e5, 1.0e5, 4.0e4, 1.0e4, 5.5e3, 3.0e3, 1.8e3, 1.1e3,
     6.5e2, 4.0e2, 2.4e2, 1.4e2, 7.0e1, 3.4e1, 1.5e1, 6.0e0, 1.8e0, 3.0e-1]
)
_LOG_ALPHA_INTERP = np.log10(_ALPHA_ANCHORS_CM1)


def alpha_silicon(wavelength_nm):
    """
    Approximate absorption coefficient of crystalline silicon, alpha(lambda),
    in 1/cm, obtained by log-linear interpolation of literature anchor
    points. Valid (and only intended for teaching use) in 250-1200 nm.
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    log_alpha = np.interp(
        wavelength_nm, _ALPHA_ANCHORS_NM, _LOG_ALPHA_INTERP,
        left=_LOG_ALPHA_INTERP[0], right=_LOG_ALPHA_INTERP[-1]
    )
    return 10 ** log_alpha


def absorption_length_um(wavelength_nm):
    """Absorption length L_alpha = 1/alpha, in micrometres."""
    alpha_cm1 = alpha_silicon(wavelength_nm)
    return 1.0 / alpha_cm1 * 1e4  # cm -> um


# ---------------------------------------------------------------------------
# Simplified AM1.5G-like reference spectrum
# ---------------------------------------------------------------------------
def am15g_simplified(wavelength_nm):
    """
    Smoothed analytic approximation of the AM1.5G solar spectral irradiance,
    in W / (m^2 nm). This is a teaching approximation (a scaled ~5778 K
    blackbody shape with the UV/near-IR roll-off of the terrestrial
    spectrum), NOT the tabulated IEC 60904-3 reference spectrum. Do not use
    for real Jsc/EQE calibration, use tabulated AM1.5G data instead.
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    wl_m = wavelength_nm * 1e-9
    T_sun = 5778.0  # K, effective blackbody temperature

    # Planck spectral radiance (per unit wavelength), arbitrary units
    with np.errstate(over="ignore", divide="ignore"):
        num = 2 * H * C0 ** 2 / wl_m ** 5
        denom = np.expm1(H * C0 / (wl_m * KB * T_sun))
        spectral_radiance = num / denom

    # Empirical terrestrial roll-off: suppresses UV (< ~320 nm, ozone) and
    # softens the far edge, loosely mimicking atmospheric attenuation.
    uv_cutoff = 1 / (1 + np.exp(-(wavelength_nm - 320) / 25))
    ir_taper = np.exp(-((np.maximum(wavelength_nm - 1600, 0)) / 900) ** 2)
    shaped = spectral_radiance * uv_cutoff * ir_taper

    # Normalise so that the integral over 280-4000 nm matches the standard
    # AM1.5G total irradiance of ~1000 W/m^2.
    wl_grid = np.linspace(280, 4000, 4000)
    wl_grid_m = wl_grid * 1e-9
    with np.errstate(over="ignore", divide="ignore"):
        num_g = 2 * H * C0 ** 2 / wl_grid_m ** 5
        denom_g = np.expm1(H * C0 / (wl_grid_m * KB * T_sun))
        rad_g = num_g / denom_g
    uv_g = 1 / (1 + np.exp(-(wl_grid - 320) / 25))
    ir_g = np.exp(-((np.maximum(wl_grid - 1600, 0)) / 900) ** 2)
    shaped_g = rad_g * uv_g * ir_g
    total = _trapz(shaped_g, wl_grid)  # W/m^2 (per the arbitrary radiance units)
    scale = 1000.0 / total

    return shaped * scale  # W / (m^2 nm)


def photon_flux_spectral(wavelength_nm, irradiance_W_m2_nm):
    """
    Convert spectral irradiance E(lambda) [W/(m^2 nm)] to spectral photon
    flux Phi0(lambda) [photons / (s m^2 nm)] via Phi0 = E * lambda / (h c).
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    wl_m = wavelength_nm * 1e-9
    return irradiance_W_m2_nm * wl_m / (H * C0)


# ---------------------------------------------------------------------------
# Carrier generation (Lambert-Beer)
# ---------------------------------------------------------------------------
def generation_profile(z_um, wavelength_nm, R=0.0):
    """
    Normalised carrier generation rate per unit depth,
    g(z, lambda) = (1 - R) * alpha * exp(-alpha z), in 1/um, such that
    integrating g(z) dz over 0..infinity gives (1 - R).

    z_um : depth into the cell, in micrometres (0 = illuminated surface).
    wavelength_nm : photon wavelength, in nm.
    R : reflectance at this wavelength (0-1).
    """
    z_um = np.asarray(z_um, dtype=float)
    alpha_cm1 = alpha_silicon(wavelength_nm)
    alpha_um1 = alpha_cm1 * 1e-4  # 1/cm -> 1/um
    return (1 - R) * alpha_um1 * np.exp(-alpha_um1 * z_um)


# ---------------------------------------------------------------------------
# Collection efficiency (base region), Eq. (10)-(11) of the lab manual
# ---------------------------------------------------------------------------
def effective_diffusion_length_um(L_um, W_um, S_cm_s, D_cm2_s):
    """
    Effective diffusion length L_eff, accounting for rear surface
    recombination velocity S (Eq. 11):

        L_eff = L * [S sinh(W/L) + D cosh(W/L)] / [S cosh(W/L) + D sinh(W/L)]

    L_um, W_um : diffusion length and cell/base thickness, in micrometres.
    S_cm_s : rear surface recombination velocity, cm/s.
    D_cm2_s : minority-carrier diffusion constant, cm^2/s.
    """
    L_cm = L_um * 1e-4
    W_cm = W_um * 1e-4
    x = W_cm / L_cm
    num = S_cm_s * np.sinh(x) + D_cm2_s * np.cosh(x)
    den = S_cm_s * np.cosh(x) + D_cm2_s * np.sinh(x)
    L_eff_cm = L_cm * num / den
    return L_eff_cm * 1e4  # back to um


def collection_efficiency(z_um, L_um, W_um, S_cm_s, D_cm2_s=27.0):
    """
    Collection probability eta_c(z) in the base region (Eq. 10):

        eta_c(z) = cosh(z/L) - (L/L_eff) * sinh(z/L)

    z_um : depth, micrometres (0 = junction / start of base region).
    L_um : (bulk) minority-carrier diffusion length, micrometres.
    W_um : base thickness, micrometres.
    S_cm_s : rear surface recombination velocity, cm/s.
    D_cm2_s : minority-carrier diffusion constant, cm^2/s (default: electrons
              in p-type Si, ~27 cm^2/s).
    """
    z_um = np.asarray(z_um, dtype=float)
    L_eff_um = effective_diffusion_length_um(L_um, W_um, S_cm_s, D_cm2_s)
    return np.cosh(z_um / L_um) - (L_um / L_eff_um) * np.sinh(z_um / L_um)


def front_surface_transmission(wavelength_nm, edge_nm=380.0, width_nm=45.0):
    """
    Phenomenological factor (0-1) for the fraction of light entering the
    cell that is NOT lost to parasitic absorption and recombination in the
    front "dead layer" (ARC + heavily doped emitter). Real short-wavelength
    photons are absorbed within tens of nm of the surface (see
    `absorption_length_um`), in a region with poor carrier collection; the
    base-region `collection_efficiency` model above assumes ideal (100%)
    collection right at its own z=0 (the front junction), so this factor
    supplies the missing front-region loss as a smooth step, rising from
    ~0 at deep UV to ~1 by the visible range. This reproduces the well
    known suppressed "blue response" of real EQE curves (cf. Fig. 5 of the
    lab manual) without modelling the emitter transport explicitly.
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    return 1.0 / (1.0 + np.exp(-(wavelength_nm - edge_nm) / (width_nm / 4.0)))


# ---------------------------------------------------------------------------
# EQE spectrum from the 1-D generation/collection model
# ---------------------------------------------------------------------------
def eqe_spectrum(wavelength_nm, W_um=160.0, L_um=200.0, S_cm_s=100.0,
                  D_cm2_s=27.0, reflectance=None, front_edge_nm=380.0,
                  front_width_nm=45.0, n_z=400):
    """
    External quantum efficiency EQE(lambda) obtained by integrating the
    generation-collection model over the base thickness, and applying a
    front dead-layer loss factor at short wavelengths:

        EQE(lambda) = T_front(lambda) * (1 - R(lambda)) *
                       integral_0^W  alpha(lambda) * exp(-alpha(lambda) z)
                       * eta_c(z) dz

    The depth integral is evaluated with adaptive quadrature (not a fixed
    grid): at short wavelengths the absorption length can be tens of nm,
    far finer than a uniform grid spanning a ~160 um cell would resolve,
    which would otherwise silently overestimate the integral.

    wavelength_nm : array of wavelengths, nm.
    W_um : base thickness, micrometres (typical Si wafer cell: ~160 um).
    L_um : bulk diffusion length, micrometres.
    S_cm_s : rear surface recombination velocity, cm/s.
    D_cm2_s : diffusion constant, cm^2/s.
    reflectance : None, a scalar, or an array matching wavelength_nm; if
                  None, a simple ARC reflectance model is used (see
                  `arc_reflectance`).
    front_edge_nm, front_width_nm : parameters of the front dead-layer loss
                  factor, see `front_surface_transmission`. Set
                  front_edge_nm=0 to disable this loss (ideal front surface).
    n_z : unused, kept for backward-compatible call signatures.
    """
    wavelength_nm = np.atleast_1d(np.asarray(wavelength_nm, dtype=float))
    if reflectance is None:
        R = arc_reflectance(wavelength_nm)
    else:
        R = np.broadcast_to(np.asarray(reflectance, dtype=float), wavelength_nm.shape)

    T_front = front_surface_transmission(wavelength_nm, edge_nm=front_edge_nm,
                                          width_nm=front_width_nm)

    eqe = np.zeros_like(wavelength_nm)

    for i, wl in enumerate(wavelength_nm):
        alpha_um1 = alpha_silicon(wl) * 1e-4  # 1/cm -> 1/um

        def integrand(z, wl=wl, R_i=R[i], alpha_um1=alpha_um1):
            g = (1 - R_i) * alpha_um1 * np.exp(-alpha_um1 * z)
            eta_c = collection_efficiency(z, L_um, W_um, S_cm_s, D_cm2_s)
            return g * eta_c

        # Split the integration range at a few absorption lengths so quad
        # resolves the (possibly very sharp) exponential near z=0 even when
        # W_um is orders of magnitude larger than 1/alpha.
        breakpoint_um = min(W_um, 10.0 / alpha_um1) if alpha_um1 > 0 else W_um
        val1, _ = quad(integrand, 0, breakpoint_um, limit=200)
        val2 = 0.0
        if breakpoint_um < W_um:
            val2, _ = quad(integrand, breakpoint_um, W_um, limit=200)
        eqe[i] = T_front[i] * (val1 + val2)

    return np.clip(eqe, 0.0, 1.0)


def arc_reflectance(wavelength_nm, R_min=0.03, lambda_min_nm=600.0,
                     width_nm=250.0, R_uv=0.35, R_rear_onset_nm=950.0,
                     R_rear_max=0.12):
    """
    Simple reflectance model for a solar cell with an anti-reflection
    coating (ARC): a Gaussian-shaped dip centred at `lambda_min_nm` (the
    ARC design wavelength) rising towards `R_uv` in the UV, plus a gradual
    increase above `R_rear_onset_nm` representing weak absorption / escaped
    light reaching the rear surface at long wavelengths (cf. Fig. 5 of the
    lab manual). This is a shape model for teaching, not a thin-film optics
    calculation.
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    dip = R_uv - (R_uv - R_min) * np.exp(
        -((wavelength_nm - lambda_min_nm) / width_nm) ** 2
    )
    rear_rise = R_rear_max * np.clip(
        (wavelength_nm - R_rear_onset_nm) / (LAMBDA_G_NM - R_rear_onset_nm), 0, 1
    ) ** 2
    return np.clip(dip + rear_rise, 0, 1)


# ---------------------------------------------------------------------------
# Short-circuit current density from EQE (Eq. 16-17)
# ---------------------------------------------------------------------------
def jsc_from_eqe(wavelength_nm, eqe, spectrum_W_m2_nm=None):
    """
    Short-circuit current density Jsc = q * integral Phi0(lambda) EQE(lambda) dlambda,
    in mA/cm^2.

    wavelength_nm, eqe : matching arrays.
    spectrum_W_m2_nm : spectral irradiance array [W/(m^2 nm)] matching
                        wavelength_nm; defaults to `am15g_simplified`.
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    eqe = np.asarray(eqe, dtype=float)
    if spectrum_W_m2_nm is None:
        spectrum_W_m2_nm = am15g_simplified(wavelength_nm)
    phi0 = photon_flux_spectral(wavelength_nm, spectrum_W_m2_nm)  # photons/(s m^2 nm)
    integrand = phi0 * eqe  # photons/(s m^2 nm)
    jsc_A_m2 = Q * _trapz(integrand, wavelength_nm)  # A/m^2
    return jsc_A_m2 * 1e-1  # A/m^2 -> mA/cm^2  (1 A/m^2 = 0.1 mA/cm^2)


# ---------------------------------------------------------------------------
# Spectral responsivity <-> EQE conversion (Eq. 27-28)
# ---------------------------------------------------------------------------
def responsivity_from_eqe(wavelength_nm, eqe):
    """Spectral responsivity s(lambda) = EQE * q * lambda / (h c), in A/W."""
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    wl_m = wavelength_nm * 1e-9
    return np.asarray(eqe, dtype=float) * Q * wl_m / (H * C0)


def eqe_from_responsivity(wavelength_nm, responsivity_A_W):
    """Inverse of `responsivity_from_eqe`: EQE(lambda) = s(lambda) * h c / (q lambda)."""
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    wl_m = wavelength_nm * 1e-9
    return np.asarray(responsivity_A_W, dtype=float) * H * C0 / (Q * wl_m)


# ---------------------------------------------------------------------------
# Synthetic "measurement" data generators, for illustrating the
# differential-spectral-responsivity (DSR) measurement workflow.
# ---------------------------------------------------------------------------
def synthetic_dsr_measurement(wavelength_nm, eqe_true, noise_level=0.01,
                               calibration_scale=1.0, random_state=None):
    """
    Generate a synthetic "measured" differential spectral responsivity
    curve from a true EQE spectrum: converts to responsivity, applies a
    (relative) calibration scaling factor C_ref/C_test (see Eq. 22-25 of
    the lab manual) and adds Gaussian noise.

    Returns the noisy, mis-scaled responsivity array s_tilde(lambda), in A/W.
    """
    rng = np.random.default_rng(random_state)
    s_true = responsivity_from_eqe(wavelength_nm, eqe_true)
    noise = rng.normal(0, noise_level * np.max(s_true), size=s_true.shape)
    return calibration_scale * s_true + noise


def bias_ramp_dsr(E_bias, s_stc, nonlinearity=0.0):
    """
    Model the differential spectral responsivity s_tilde as a function of
    bias irradiance E_bias (W/m^2), at fixed wavelength.

    For an ideal, linear solar cell, s_tilde is independent of E_bias
    (`nonlinearity=0`). A nonlinear cell shows a systematic increase or
    decrease of s_tilde with bias level, saturating at high irradiance,
    modelled here as:

        s_tilde(E_bias) = s_stc * (1 + nonlinearity * (1 - exp(-E_bias / 300)))

    E_bias : bias irradiance, W/m^2 (typically evaluated 0-1000 W/m^2).
    s_stc  : the true DSR under standard test conditions (E_bias -> 1000 W/m^2), A/W.
    nonlinearity : relative deviation at high bias irradiance (e.g. 0.15 = 15%).
    """
    E_bias = np.asarray(E_bias, dtype=float)
    return s_stc * (1 + nonlinearity * (1 - np.exp(-E_bias / 300.0)))
