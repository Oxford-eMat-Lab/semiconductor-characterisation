"""
tlm_helper.py
-------------
Helper functions for the Transfer Length Method (TLM) teaching notebook,
`tlm_analysis.ipynb`. Equation numbers quoted in the docstrings below,
e.g. "Eq. (5)", refer to the numbered equations in that notebook.

Scope
-----
These functions implement the standard 1-D transmission-line description of
a planar metal/semiconductor contact: a uniform sheet of resistance R_S
underneath contacts that inject current through a uniform specific contact
resistivity rho_C. That model is what the TLM extraction assumes, so using
it here makes the assumptions - and the ways they fail - explicit.

They are written for teaching, not for precision metrology:
  - the semiconductor layer is assumed uniform and purely ohmic,
  - the contacts are assumed identical and infinitely conductive,
  - current flow is 1-D along the strip; edge/lateral spreading outside the
    contact width is ignored (no end-effect correction),
  - the synthetic "measurements" add Gaussian noise only; real datasets also
    carry probe placement error, self-heating and non-ohmic behaviour.

For quantitative work, use a full TLM test structure with a measured
contact width, correct for probe/lead resistance, and check linearity of
every I-V curve before fitting.

UNIT CONVENTION
---------------
Lengths that describe the *test structure* (contact spacing L, contact
width W, contact length d) are in millimetres, because that is what is
drawn on a mask and what a probe station reports. Material quantities keep
their conventional units: R_S in ohm/square, rho_C in ohm.cm^2, rho in
ohm.cm, thickness in cm. Conversions happen inside these functions, never
in the notebook, so the notebook never juggles units.
"""

import re
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Unit conversions (kept explicit rather than sprinkled through the code)
# ---------------------------------------------------------------------------
MM_PER_CM = 10.0        # 1 cm = 10 mm
CM_PER_MM = 0.1         # 1 mm = 0.1 cm


def mm_to_cm(x_mm):
    """Millimetres -> centimetres."""
    return np.asarray(x_mm, dtype=float) * CM_PER_MM


def cm_to_mm(x_cm):
    """Centimetres -> millimetres."""
    return np.asarray(x_cm, dtype=float) * MM_PER_CM


# ---------------------------------------------------------------------------
# 1. Sheet resistance of a thin conducting layer
# ---------------------------------------------------------------------------
def sheet_resistance(rho_ohm_cm, thickness_cm):
    """
    Sheet resistance R_S = rho / t, in ohm/square - Eq. (2).

    rho_ohm_cm : bulk resistivity of the layer, ohm.cm.
    thickness_cm : layer thickness, cm (note: 100 nm = 100e-7 cm).

    A "square" is dimensionless: any square patch of the layer has the same
    resistance between two opposite edges, whatever its size, which is why
    R_S describes the layer independently of how it is patterned.
    """
    return np.asarray(rho_ohm_cm, dtype=float) / np.asarray(thickness_cm, dtype=float)


def strip_resistance(Rs_ohm_sq, L_mm, W_mm):
    """
    Resistance of a rectangular patch of sheet, R = R_S * L / W - Eq. (3).

    Rs_ohm_sq : sheet resistance, ohm/square.
    L_mm : length of the current path (here: the contact spacing), mm.
    W_mm : width of the strip / contact, mm.

    Only the *ratio* L/W (the "number of squares") matters, so the units of
    L and W cancel as long as they are the same - mm here.
    """
    return np.asarray(Rs_ohm_sq, dtype=float) * np.asarray(L_mm, dtype=float) / float(W_mm)


# ---------------------------------------------------------------------------
# 2. Transfer length and specific contact resistivity
# ---------------------------------------------------------------------------
def transfer_length_mm(rho_c_ohm_cm2, Rs_ohm_sq):
    """
    Transfer length L_T = sqrt(rho_C / R_S), returned in mm - Eq. (7).

    rho_c_ohm_cm2 : specific contact resistivity, ohm.cm^2.
    Rs_ohm_sq : sheet resistance, ohm/square.

    sqrt(ohm.cm^2 / (ohm/square)) = cm, which is then converted to mm.
    Physically L_T is the average distance a carrier travels laterally in
    the semiconductor *under* the contact before crossing into the metal:
    it is set by the competition between lateral transport (R_S) and
    vertical transport across the interface (rho_C).
    """
    LT_cm = np.sqrt(np.asarray(rho_c_ohm_cm2, dtype=float) /
                    np.asarray(Rs_ohm_sq, dtype=float))
    return cm_to_mm(LT_cm)


def contact_resistivity(LT_mm, Rs_ohm_sq):
    """
    Specific contact resistivity rho_C = R_S * L_T^2, in ohm.cm^2 - Eq. (7)
    solved for rho_C, using the definition rho_C = R_C * A_eff of Eq. (5).

    This is the form to prefer when reporting a TLM result, because it uses
    the two quantities the fit actually determines (slope -> R_S, and the
    intercept ratio -> L_T) rather than R_C, whose value depends on the
    contact size.
    """
    LT_cm = mm_to_cm(LT_mm)
    return np.asarray(Rs_ohm_sq, dtype=float) * LT_cm ** 2


def contact_resistance(Rs_ohm_sq, LT_mm, W_mm, d_mm=None):
    """
    Resistance of a single contact, in ohm - Eq. (9).

    Long-contact limit (d >> L_T, the usual TLM design):

        R_C = rho_C / (L_T W) = R_S L_T / W

    If the contact length `d_mm` (the dimension along the current
    direction) is given, the full transmission-line result is returned
    instead:

        R_C = (R_S L_T / W) * coth(d / L_T)

    The coth factor is within 1% of 1 once d > ~3 L_T, which is why TLM
    patterns are drawn with generously long contacts: a short contact
    inflates R_C and makes rho_C look worse than it is.
    """
    Rs_ohm_sq = np.asarray(Rs_ohm_sq, dtype=float)
    LT_mm = np.asarray(LT_mm, dtype=float)
    Rc = Rs_ohm_sq * LT_mm / float(W_mm)
    if d_mm is not None:
        Rc = Rc / np.tanh(np.asarray(d_mm, dtype=float) / LT_mm)
    return Rc


def total_resistance(L_mm, Rs_ohm_sq, LT_mm, W_mm):
    """
    The TLM master equation - Eq. (10):

        R_T(L) = (R_S / W) * (L + 2 L_T)

    i.e. a straight line in the contact spacing L with
        slope       = R_S / W
        y-intercept = 2 R_C  = 2 R_S L_T / W
        x-intercept = -2 L_T

    Those three readings off one plot are the whole method.
    """
    return (np.asarray(Rs_ohm_sq, dtype=float) / float(W_mm)) * \
           (np.asarray(L_mm, dtype=float) + 2.0 * np.asarray(LT_mm, dtype=float))


def current_under_contact(x_mm, LT_mm, I0=1.0):
    """
    Current still flowing in the semiconductor a distance x into the
    contact, I(x) = I0 exp(-x / L_T) - Eq. (6).

    x is measured from the contact edge where the current enters. This is
    the exponential decay that makes "contact area" ambiguous in a planar
    geometry: current does not enter uniformly over the whole contact
    footprint, it crowds within roughly one L_T of the leading edge.
    """
    return float(I0) * np.exp(-np.asarray(x_mm, dtype=float) /
                              np.asarray(LT_mm, dtype=float))


def effective_contact_area_cm2(LT_mm, W_mm):
    """
    Effective current-transfer area A_eff = L_T * W, in cm^2 - Eq. (8).

    The area that enters rho_C = R_C * A_eff, Eq. (5). Note it uses L_T, not the
    drawn contact length: enlarging a contact beyond a few L_T adds
    footprint but no extra current path.
    """
    return mm_to_cm(LT_mm) * mm_to_cm(W_mm)


# ---------------------------------------------------------------------------
# 3. Forward model: I-V curves of a TLM pair
# ---------------------------------------------------------------------------
def iv_curve(V, Rs_ohm_sq, L_mm, W_mm, LT_mm=0.0, R_series_ohm=0.0):
    """
    Current through one contact pair of a TLM structure, in A - Eq. (10).

    An ohmic device, so I = V / R_T with R_T from `total_resistance`, plus
    an optional `R_series_ohm` standing in for probe, lead and metal
    resistance (the R_m term dropped in Eq. (1)).

    V : applied voltage, V (scalar or array).
    LT_mm = 0 reduces this to a contact-free strip, which is how the
    notebook shows what R_S alone does before contacts are introduced.
    """
    R = total_resistance(L_mm, Rs_ohm_sq, LT_mm, W_mm) + float(R_series_ohm)
    return np.asarray(V, dtype=float) / R


def iv_curve_schottky(V, R_ohm, I_sat=1e-6, n=2.0, T=300.0):
    """
    I-V curve of a *non-ohmic* (rectifying) contact pair, for the section
    on data quality. Two back-to-back diodes in series with the sheet
    resistance give a symmetric but S-shaped curve, whose "resistance"
    depends on the voltage range you happen to fit.

    Not part of the TLM model - included precisely because fitting a
    straight line to this is the single most common way a TLM extraction
    goes silently wrong.
    """
    kT_q = 1.380649e-23 * float(T) / 1.602176634e-19
    V = np.asarray(V, dtype=float)
    # Symmetric back-to-back diode: current limited by the reverse-biased one.
    I_diode = np.sign(V) * I_sat * (np.exp(np.abs(V) / (2 * n * kT_q)) - 1.0)
    # Series sheet resistance limits the current at high bias.
    return np.sign(V) * np.minimum(np.abs(I_diode), np.abs(V) / R_ohm)


# ---------------------------------------------------------------------------
# 4. Synthetic measurements
# ---------------------------------------------------------------------------
def synthetic_iv_measurement(V, R_true_ohm, noise_pA=0.0, offset_uV=0.0,
                             seed=None):
    """
    One measured I-V sweep of a resistor: Ohm's law plus Gaussian current
    noise and an optional voltage offset (a real source-meter has both).

    Returns the current array, A.
    """
    rng = np.random.default_rng(seed)
    V = np.asarray(V, dtype=float)
    I = (V + offset_uV * 1e-6) / float(R_true_ohm)
    return I + rng.normal(0.0, noise_pA * 1e-12, size=V.shape)


def synthetic_tlm_dataset(spacings_mm, Rs_ohm_sq=250.0, rho_c_ohm_cm2=1e-3,
                          W_mm=10.0, n_repeats=4, V_max=0.2, n_points=101,
                          noise_pA=200.0, spacing_error_um=8.0,
                          R_series_ohm=0.0, seed=0):
    """
    A complete synthetic TLM measurement: for each contact spacing, several
    nominally identical structures are "measured" over a voltage sweep.

    The two error sources are deliberately different in kind:
      - `noise_pA`   : random current noise, which averages away with
                       repeats and mostly affects each I-V fit;
      - `spacing_error_um` : a lithography/probe-placement error on the
                       actual spacing, which does *not* average away and is
                       what usually dominates the intercept - and therefore
                       rho_C.

    Returns a dict with:
      'V'          : voltage array, V, shape (n_points,)
      'I'          : currents, A, shape (n_spacings, n_repeats, n_points)
      'spacings_mm': the nominal (drawn) spacings, mm
      'truth'      : dict of the R_S, rho_C, L_T, R_C used to generate it,
                     so the notebook can check what the extraction recovers.
    """
    rng = np.random.default_rng(seed)
    spacings_mm = np.asarray(spacings_mm, dtype=float)
    V = np.linspace(-V_max, V_max, n_points)

    LT_mm = transfer_length_mm(rho_c_ohm_cm2, Rs_ohm_sq)
    Rc = contact_resistance(Rs_ohm_sq, LT_mm, W_mm)

    I = np.empty((spacings_mm.size, int(n_repeats), V.size))
    for i, L in enumerate(spacings_mm):
        for j in range(int(n_repeats)):
            L_actual = L + rng.normal(0.0, spacing_error_um * 1e-3)
            R_true = total_resistance(L_actual, Rs_ohm_sq, LT_mm, W_mm) + R_series_ohm
            I[i, j] = V / R_true + rng.normal(0.0, noise_pA * 1e-12, size=V.size)

    return {
        'V': V,
        'I': I,
        'spacings_mm': spacings_mm,
        'W_mm': W_mm,
        'truth': {
            'Rs_ohm_sq': Rs_ohm_sq,
            'rho_c_ohm_cm2': rho_c_ohm_cm2,
            'LT_mm': float(LT_mm),
            'Rc_ohm': float(Rc),
        },
    }


# ---------------------------------------------------------------------------
# 5. Step 1 of the extraction: I-V sweep -> one resistance
# ---------------------------------------------------------------------------
def resistance_from_iv(V, I, V_window=None):
    """
    Resistance from the slope of a measured I-V sweep - Eq. (11).

    Fits I = aV + b and returns 1/a. Fitting the *slope* rather than
    dividing V by I point-by-point removes any current offset b (leakage,
    thermal EMF) from the answer.

    V_window : optional (Vmin, Vmax); restrict the fit to a symmetric
        window around zero bias. Widening this window on a non-ohmic
        device changes the answer - which is exactly the diagnostic used
        in the data-quality section.
    """
    V = np.asarray(V, dtype=float)
    I = np.asarray(I, dtype=float)
    good = np.isfinite(V) & np.isfinite(I)
    if V_window is not None:
        good &= (V >= V_window[0]) & (V <= V_window[1])
    if good.sum() < 2:
        return np.nan
    a, _b = np.polyfit(V[good], I[good], 1)
    return 1.0 / a if a != 0 else np.nan


def iv_linearity(V, I, V_window=None):
    """
    Coefficient of determination R^2 of the straight-line fit to an I-V
    sweep - a one-number check that a contact pair is ohmic before its
    resistance is used in a TLM fit.

    R^2 below ~0.999 on a clean sweep means the curve is bending; treat the
    extracted resistance with suspicion.
    """
    V = np.asarray(V, dtype=float)
    I = np.asarray(I, dtype=float)
    good = np.isfinite(V) & np.isfinite(I)
    if V_window is not None:
        good &= (V >= V_window[0]) & (V <= V_window[1])
    if good.sum() < 3:
        return np.nan
    V, I = V[good], I[good]
    a, b = np.polyfit(V, I, 1)
    resid = I - (a * V + b)
    ss_tot = np.sum((I - I.mean()) ** 2)
    return 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan


def resistances_per_spacing(dataset, V_window=None):
    """
    Reduce a `synthetic_tlm_dataset` (or a dict of the same shape built by
    `load_tlm_folder`) to one mean resistance and one standard deviation
    per contact spacing.

    Returns (spacings_mm, R_mean_ohm, R_std_ohm, R_all_ohm) where R_all has
    shape (n_spacings, n_repeats) so individual structures can still be
    inspected or rejected.
    """
    V = dataset['V']
    I = dataset['I']
    n_s, n_r = I.shape[0], I.shape[1]
    R_all = np.empty((n_s, n_r))
    for i in range(n_s):
        for j in range(n_r):
            R_all[i, j] = resistance_from_iv(V, I[i, j], V_window=V_window)
    R_mean = np.nanmean(R_all, axis=1)
    R_std = np.nanstd(R_all, axis=1, ddof=1) if n_r > 1 else np.zeros(n_s)
    return dataset['spacings_mm'], R_mean, R_std, R_all


# ---------------------------------------------------------------------------
# 6. Step 2 of the extraction: R_T vs L -> R_S, R_C, L_T, rho_C
# ---------------------------------------------------------------------------
class TLMResult(dict):
    """
    Result of a TLM fit. A plain dict (so it prints and unpacks predictably)
    with attribute access for readability: `res.rho_c_ohm_cm2` works, and so
    does `res['rho_c_ohm_cm2']`.

    Keys: slope_ohm_per_mm, intercept_ohm, Rs_ohm_sq, Rc_ohm, LT_mm,
    rho_c_ohm_cm2, and the matching *_err one-sigma uncertainties, plus
    r_squared, L_mm, R_ohm, W_mm.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def summary(self):
        """A formatted multi-line report of the extracted parameters."""
        lines = [
            "TLM extraction",
            "--------------",
            f"  fit           R_T = {self['slope_ohm_per_mm']:.4g} * L + "
            f"{self['intercept_ohm']:.4g}   (R^2 = {self['r_squared']:.5f})",
            f"  contact width W       = {self['W_mm']:.3g} mm",
            f"  sheet resist. R_S     = {self['Rs_ohm_sq']:.4g} "
            f"+/- {self['Rs_err']:.2g} ohm/sq",
            f"  contact res.  R_C     = {self['Rc_ohm']:.4g} "
            f"+/- {self['Rc_err']:.2g} ohm",
            f"  transfer len. L_T     = {self['LT_mm']:.4g} "
            f"+/- {self['LT_err']:.2g} mm",
            f"  contact resy. rho_C   = {self['rho_c_ohm_cm2']:.4g} "
            f"+/- {self['rho_c_err']:.2g} ohm.cm^2",
        ]
        return "\n".join(lines)

    def __str__(self):
        return self.summary()


def fit_tlm(L_mm, R_ohm, W_mm, sigma_ohm=None):
    """
    The TLM extraction itself - Eqs. (12)-(15).

    Fits R_T = m L + b to the measured (spacing, resistance) pairs, then:

        R_S   = m * W                      Eq. (12)
        R_C   = b / 2                      Eq. (13)
        L_T   = b / (2 m)                  Eq. (14)
        rho_C = R_S * L_T^2                Eq. (15)

    Uncertainties are propagated from the covariance matrix of the linear
    least-squares fit, including the m-b correlation, which matters: m and
    b from a fit over a limited range of L are strongly anti-correlated, so
    treating their errors as independent overstates the error on L_T.

    L_mm : contact spacings, mm.
    R_ohm : total measured resistance at each spacing, ohm.
    W_mm : contact width, mm (measure it; a wrong W scales R_S and rho_C).
    sigma_ohm : optional per-point one-sigma uncertainties, used to weight
        the fit. Pass the spread over repeated structures.

    Returns a `TLMResult`.
    """
    L = np.asarray(L_mm, dtype=float)
    R = np.asarray(R_ohm, dtype=float)
    good = np.isfinite(L) & np.isfinite(R)
    L, R = L[good], R[good]
    if L.size < 3:
        raise ValueError("need at least 3 spacings for a meaningful TLM fit "
                         f"(got {L.size})")

    w = None
    if sigma_ohm is not None:
        s = np.asarray(sigma_ohm, dtype=float)[good]
        if np.all(np.isfinite(s)) and np.all(s > 0):
            w = 1.0 / s

    (m, b), cov = np.polyfit(L, R, 1, w=w, cov=True)
    var_m, var_b, cov_mb = cov[0, 0], cov[1, 1], cov[0, 1]
    err_m, err_b = np.sqrt(var_m), np.sqrt(var_b)

    resid = R - (m * L + b)
    ss_tot = np.sum((R - R.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan

    Rs = m * W_mm
    Rs_err = err_m * W_mm

    Rc = b / 2.0
    Rc_err = err_b / 2.0

    LT = b / (2.0 * m)
    # d(LT)/db = 1/(2m), d(LT)/dm = -b/(2m^2), including the covariance term
    dLT_db = 1.0 / (2.0 * m)
    dLT_dm = -b / (2.0 * m ** 2)
    LT_var = (dLT_db ** 2 * var_b + dLT_dm ** 2 * var_m
              + 2.0 * dLT_db * dLT_dm * cov_mb)
    LT_err = np.sqrt(max(LT_var, 0.0))

    rho_c = contact_resistivity(LT, Rs)
    # rho_C = Rs * LT^2 -> relative errors add as (dRs/Rs) and 2*(dLT/LT)
    rel = np.sqrt((Rs_err / Rs) ** 2 + (2.0 * LT_err / LT) ** 2)
    rho_c_err = abs(rho_c) * rel

    return TLMResult(
        slope_ohm_per_mm=float(m), slope_err=float(err_m),
        intercept_ohm=float(b), intercept_err=float(err_b),
        Rs_ohm_sq=float(Rs), Rs_err=float(Rs_err),
        Rc_ohm=float(Rc), Rc_err=float(Rc_err),
        LT_mm=float(LT), LT_err=float(LT_err),
        rho_c_ohm_cm2=float(rho_c), rho_c_err=float(rho_c_err),
        r_squared=float(r2), W_mm=float(W_mm),
        L_mm=L, R_ohm=R,
    )


def tlm_validity_report(result, warn_ratio=0.2, min_r2=0.995):
    """
    Check a `TLMResult` against the conditions the method assumes, and
    return a list of human-readable warnings (empty list = nothing flagged).

    Checks, in order of how often they bite:
      1. L_T must be small compared with the largest measured spacing.
         The intercept is an *extrapolation* to L = 0; if L_T approaches
         L_max the extrapolation dominates the answer.
      2. The intercept must be positive. A negative one is unphysical
         (it implies negative contact resistance) and means either a
         non-linear dataset or a systematic spacing error.
      3. The fit must actually be linear (R^2).
      4. At least 4 spacings, spanning a decent range, for the fit to be
         constrained rather than merely satisfied.
    """
    msgs = []
    L = np.asarray(result['L_mm'], dtype=float)
    L_max = L.max()
    LT = result['LT_mm']

    if not np.isfinite(LT):
        msgs.append("L_T is not finite - the fit failed.")
    elif LT <= 0:
        msgs.append(
            f"Negative transfer length (L_T = {LT:.3g} mm): the fitted "
            "intercept is below zero, which is unphysical. Suspect a "
            "non-ohmic contact, a spacing offset, or an uncorrected series "
            "resistance.")
    elif LT > warn_ratio * L_max:
        msgs.append(
            f"L_T = {LT:.3g} mm is {LT / L_max:.0%} of the largest spacing "
            f"({L_max:.3g} mm). The intercept is then mostly extrapolation; "
            "rho_C from this fit is unreliable. Measure larger spacings, so "
            "that L_max is several times L_T.")

    if result['r_squared'] < min_r2:
        msgs.append(
            f"R^2 = {result['r_squared']:.4f} < {min_r2}: R_T vs L is not "
            "straight. Check every I-V sweep for non-ohmic behaviour before "
            "trusting the intercept.")

    if L.size < 4:
        msgs.append(
            f"Only {L.size} spacings. Four or more, spanning at least a "
            "factor of 5 in L, make the intercept far better constrained.")
    elif L_max / L.min() < 3:
        msgs.append(
            f"Spacings span only a factor {L_max / L.min():.1f}. A short "
            "lever arm makes the extrapolated intercept noisy.")

    return msgs


def monte_carlo_uncertainty(L_mm, R_ohm, sigma_ohm, W_mm, n_trials=2000,
                            seed=1):
    """
    Repeat the TLM fit on `n_trials` noisy replicas of the same dataset and
    return the resulting distributions.

    Why bother when `fit_tlm` already propagates errors: rho_C depends on
    the *square* of the intercept ratio, so its error distribution is
    skewed rather than Gaussian, and the linear propagation understates the
    upper tail. Seeing the histogram is more honest than a +/- number.

    Returns a dict with arrays 'Rs_ohm_sq', 'Rc_ohm', 'LT_mm',
    'rho_c_ohm_cm2', one entry per successful trial.
    """
    rng = np.random.default_rng(seed)
    L = np.asarray(L_mm, dtype=float)
    R = np.asarray(R_ohm, dtype=float)
    s = np.broadcast_to(np.asarray(sigma_ohm, dtype=float), R.shape)

    out = {'Rs_ohm_sq': [], 'Rc_ohm': [], 'LT_mm': [], 'rho_c_ohm_cm2': []}
    for _ in range(int(n_trials)):
        R_try = R + rng.normal(0.0, s)
        try:
            res = fit_tlm(L, R_try, W_mm)
        except (ValueError, np.linalg.LinAlgError):
            continue
        for k in out:
            out[k].append(res[k])
    return {k: np.asarray(v) for k, v in out.items()}


# ---------------------------------------------------------------------------
# 7. Reading a real measurement folder
# ---------------------------------------------------------------------------
def spacing_from_filename_mm(name):
    """
    Pull the contact spacing out of a file name such as
    'sample1_2.5mm_run3.csv' -> 2.5. Returns NaN if no '<number>mm' pattern
    is present, so unnamed files are skipped rather than silently mis-fit.
    """
    m = re.search(r'(\d+(?:\.\d+)?)\s*mm', str(name).lower())
    return float(m.group(1)) if m else np.nan


def load_tlm_folder(folder, W_mm, pattern='*.csv', v_col=0, i_col=1,
                    columns_per_sweep=None, skiprows=0, delimiter=','):
    """
    Load a folder of measured I-V files into the same dict layout that
    `synthetic_tlm_dataset` produces, so every downstream function in this
    module works unchanged on real data.

    Expected layout: one file per contact spacing, with the spacing in the
    file name (see `spacing_from_filename_mm`). Each file holds one or more
    voltage/current column pairs - repeated structures at that spacing.

    folder : directory to read.
    W_mm : contact width, mm.
    columns_per_sweep : column stride between successive sweeps. Defaults
        to `i_col + 1`; set it to 3 for exporters that write a blank or
        timestamp column between sweeps.
    skiprows, delimiter : passed to numpy.genfromtxt for other formats.

    Because files may hold different numbers of points, all sweeps are
    resampled onto the voltage grid of the first sweep before stacking.
    Raises FileNotFoundError if the folder holds no matching files, so a
    typo'd path fails loudly instead of returning an empty fit.
    """
    folder = Path(folder)
    files = sorted(folder.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no files matching {pattern!r} in {folder}")

    stride = columns_per_sweep or (i_col + 1)

    per_file = []
    for f in files:
        L = spacing_from_filename_mm(f.stem)
        if not np.isfinite(L):
            continue
        raw = np.genfromtxt(f, delimiter=delimiter, skip_header=skiprows or 1,
                            invalid_raise=False)
        raw = np.atleast_2d(raw)
        sweeps = []
        for c in range(0, raw.shape[1] - 1, stride):
            V = raw[:, c + v_col]
            I = raw[:, c + i_col]
            keep = np.isfinite(V) & np.isfinite(I)
            if keep.sum() >= 2:
                sweeps.append((V[keep], I[keep]))
        if sweeps:
            per_file.append((L, sweeps))

    if not per_file:
        raise FileNotFoundError(
            f"found files in {folder} but none had a '<number>mm' spacing in "
            "the name and at least one usable V/I column pair")

    per_file.sort(key=lambda t: t[0])
    V_grid = per_file[0][1][0][0]
    n_rep = max(len(s) for _, s in per_file)

    I_stack = np.full((len(per_file), n_rep, V_grid.size), np.nan)
    for i, (_L, sweeps) in enumerate(per_file):
        for j, (V, I) in enumerate(sweeps):
            order = np.argsort(V)
            I_stack[i, j] = np.interp(V_grid, V[order], I[order],
                                      left=np.nan, right=np.nan)

    return {
        'V': V_grid,
        'I': I_stack,
        'spacings_mm': np.array([L for L, _ in per_file]),
        'W_mm': float(W_mm),
        'files': [f.name for f in files],
    }


def analyse_tlm(dataset, W_mm=None, V_window=None):
    """
    End-to-end convenience wrapper: dataset -> per-spacing resistances ->
    TLM fit -> validity report.

    Returns (result, spacings_mm, R_mean, R_std, warnings). Use this once
    the individual steps are understood; the notebook walks through them
    separately first so nothing is hidden.
    """
    W = W_mm if W_mm is not None else dataset['W_mm']
    L, R_mean, R_std, _R_all = resistances_per_spacing(dataset,
                                                       V_window=V_window)
    sigma = R_std if np.all(R_std > 0) else None
    result = fit_tlm(L, R_mean, W, sigma_ohm=sigma)
    return result, L, R_mean, R_std, tlm_validity_report(result)
