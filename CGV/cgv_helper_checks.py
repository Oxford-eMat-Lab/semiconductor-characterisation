"""Analytic-limit checks for CGV/cgv_helper.py.

Not part of the repo layout (only KPSPV has a checks file so far) - run
with:  python3 cgv_helper_checks.py
"""
import numpy as np
import cgv_helper as cgv

ok = fail = 0
def check(name, cond, detail=""):
    global ok, fail
    if cond: ok += 1; print(f"  PASS  {name}   {detail}")
    else:    fail += 1; print(f"  FAIL  {name}   {detail}")

Na = 1e16  # p-type substrate, cm^-3
tox = 10.0  # nm
area = 1e-4  # cm^2 (100 um x 100 um pad)

print("== 1. oxide and Debye length ==")
Cox_pa = cgv.oxide_capacitance_per_area(tox)
check("Cox scales as 1/tox", abs(cgv.oxide_capacitance_per_area(2*tox) - Cox_pa/2) < 1e-20)
Cox_total = cgv.oxide_capacitance(area, tox)
check("Cox total = Cox_pa * area", abs(Cox_total - Cox_pa*area) < 1e-25, f"{Cox_total:.4e} F")
Ld = cgv.debye_length_cm(Na)
check("Debye length ~40 nm at 1e16 cm^-3", 3e-6 < Ld < 5e-6, f"L_D={Ld*1e7:.1f} nm")

print("\n== 2. fermi potential ==")
phiF = cgv.fermi_potential(Na, 'p')
check("phi_F positive for p-type, ~0.35 V", 0.30 < phiF < 0.40, f"phi_F={phiF:.4f} V")
phiF_n = cgv.fermi_potential(Na, 'n')
check("phi_F negative for n-type, equal magnitude", abs(phiF_n + phiF) < 1e-12)

print("\n== 3. exact Q_sc(phi_s): flat band, accumulation, depletion ==")
Q0 = cgv.space_charge_density(0.0, Na, 'p')
check("Q_sc(0) = 0 exactly (flat band)", abs(Q0) < 1e-30, f"{Q0:.3e} C/cm^2")
Qacc = cgv.space_charge_density(-0.3, Na, 'p')
check("accumulation (phi_s<0, p-type) gives POSITIVE Q_sc", Qacc > 0, f"Q_sc={Qacc:.3e} C/cm^2")
Qdep = cgv.space_charge_density(0.2, Na, 'p')
check("depletion (0<phi_s<2phiF, p-type) gives NEGATIVE Q_sc", Qdep < 0, f"Q_sc={Qdep:.3e} C/cm^2")

# depletion approximation cross-check: for small positive phi_s (well
# below inversion), the exact |Q_sc| should match the depletion-approx
# formula sqrt(2*eps_s*q*N*phi_s) to a few percent.
phi_test = 0.15
Q_exact = abs(cgv.space_charge_density(phi_test, Na, 'p'))
Q_dep_approx = np.sqrt(2*cgv.EPS_SI*cgv.Q*Na*phi_test)
rel_err = abs(Q_exact - Q_dep_approx)/Q_dep_approx
check("exact Q_sc matches depletion approx to <12% at phi_s=0.15V (phi_F=0.358V)",
      rel_err < 0.12, f"exact={Q_exact:.4e}, approx={Q_dep_approx:.4e}, rel_err={rel_err:.3%}")

print("\n== 4. semiconductor capacitance C_s(phi_s) ==")
Cs_acc = cgv.semiconductor_capacitance(-0.3, Na, 'p')
Cs_dep = cgv.semiconductor_capacitance(0.2, Na, 'p')
check("C_s much larger in accumulation than in depletion",
      Cs_acc > 20*Cs_dep, f"Cs_acc={Cs_acc:.3e}, Cs_dep={Cs_dep:.3e} F/cm^2")
Cs_fb = cgv.semiconductor_capacitance(1e-6, Na, 'p')
Cs_fb_analytic = cgv.EPS_SI/Ld
check("C_s(flat band) matches eps_s/L_D to <1%",
      abs(Cs_fb-Cs_fb_analytic)/Cs_fb_analytic < 0.01,
      f"numeric={Cs_fb:.4e}, analytic={Cs_fb_analytic:.4e} F/cm^2")

print("\n== 5. charge-balance solve: surface_potential_from_bias ==")
VFB = -0.9
phi_s_at_VFB = cgv.surface_potential_from_bias(VFB, Cox_pa, VFB, Na, 'p')
check("phi_s = 0 exactly when V_G = V_FB", abs(phi_s_at_VFB) < 1e-8, f"phi_s={float(phi_s_at_VFB):.3e} V")
phi_s_acc = cgv.surface_potential_from_bias(VFB-3.0, Cox_pa, VFB, Na, 'p')
check("large negative bias (below VFB, p-type) drives ACCUMULATION (phi_s<0)",
      phi_s_acc < 0, f"phi_s={phi_s_acc:.3f} V")
phi_s_dep = cgv.surface_potential_from_bias(VFB+3.0, Cox_pa, VFB, Na, 'p')
check("large positive bias (above VFB, p-type) drives DEPLETION/INVERSION (phi_s>0)",
      phi_s_dep > 0, f"phi_s={phi_s_dep:.3f} V")

print("\n== 6. HF/LF/deep-depletion C-V curves ==")
Vg = np.linspace(VFB-3.0, VFB+3.0, 41)
C_hf = cgv.hf_cv_curve(Vg, area, tox, Na, 'p', VFB)
C_lf = cgv.lf_cv_curve(Vg, area, tox, Na, 'p', VFB)
check("HF accumulation plateau within 3% of Cox (thin 10nm oxide doesn't fully flatten)",
      abs(C_hf[0]-Cox_total)/Cox_total < 0.03, f"C_hf(acc)={C_hf[0]:.4e}, Cox={Cox_total:.4e} F")
check("LF accumulation plateau within 3% of Cox",
      abs(C_lf[0]-Cox_total)/Cox_total < 0.03, f"C_lf(acc)={C_lf[0]:.4e} F")
check("HF inversion plateau < Cox (partial series combination)",
      C_hf[-1] < 0.9*Cox_total, f"C_hf(inv)={C_hf[-1]:.4e}, Cox={Cox_total:.4e} F")
check("LF inversion RISES BACK toward Cox (minority carriers keep up)",
      C_lf[-1] > 1.5*C_hf[-1], f"C_lf(inv)={C_lf[-1]:.4e}, C_hf(inv)={C_hf[-1]:.4e} F")
check("HF curve is monotonically non-increasing from accumulation to the HF plateau (<0.1% wiggle tolerated)",
      np.all(np.diff(C_hf) <= 1e-3*Cox_total), f"max positive step={np.diff(C_hf).max():.3e} F")

C_dd = cgv.deep_depletion_cv(Vg, area, tox, Na, 'p', VFB)
check("deep depletion falls BELOW the HF equilibrium minimum",
      C_dd[-1] < C_hf[-1], f"C_dd(end)={C_dd[-1]:.4e}, C_hf(inv)={C_hf[-1]:.4e} F")

print("\n== 7. flat-band voltage / effective oxide charge round-trip ==")
phi_ms = -0.85
Qeff_true = 2e-9  # C/cm^2
VFB_calc = cgv.flatband_voltage(phi_ms, Qeff_true, Cox_pa)
Qeff_recovered = cgv.effective_oxide_charge(phi_ms, VFB_calc, Cox_pa)
check("effective_oxide_charge inverts flatband_voltage",
      abs(Qeff_recovered - Qeff_true) < 1e-20, f"Qeff={Qeff_recovered:.4e} C/cm^2")

print("\n== 8. threshold voltage ==")
VT = cgv.threshold_voltage(Cox_pa, Na, 'p', VFB)
check("V_T > V_FB for p-type (must deplete further to invert)", VT > VFB, f"VT={VT:.3f}, VFB={VFB:.3f} V")
VT_n = cgv.threshold_voltage(Cox_pa, Na, 'n', VFB)
check("V_T < V_FB for n-type (mirror)", VT_n < VFB, f"VT_n={VT_n:.3f} V")

print("\n== 9. doping profile extraction ==")
N_true = 3e16
VFB2 = -0.8
phiF2 = cgv.fermi_potential(N_true, 'p')
phi_s_min2 = cgv._phi_s_min_capacitance(N_true, 'p')
# Stay well clear of the inversion plateau (where dC/dV -> 0 and the
# doping-profile formula is singular by construction - Sec. 9 says so).
Vg_hi = VFB2 + 0.7 * (phi_s_min2)  # comfortably inside pure depletion
Vg2 = np.linspace(VFB2+0.15, Vg_hi, 60)
C2 = cgv.hf_cv_curve(Vg2, area, tox, N_true, 'p', VFB2)
W, N_extracted = cgv.doping_profile_from_cv(Vg2, C2, area)
mid = slice(10, 45)
rel_err_N = np.abs(N_extracted[mid] - N_true) / N_true
check("uniform doping recovered within 5% in the well-resolved depletion region",
      np.median(rel_err_N) < 0.05, f"median rel. error = {np.median(rel_err_N):.3%}")

print("\n== 10. Mott-Schottky (no-oxide) reduction ==")
Vbi_true = 0.75
V_ms = np.linspace(-8, 0, 60)
C_ms = cgv.synthetic_mott_schottky(V_ms, N_true, Vbi_true, area, noise_frac=0.0)
N_fit, N_fit_err, Vbi_fit, Vbi_fit_err = cgv.mott_schottky_fit(V_ms, C_ms, area)
check("Mott-Schottky fit recovers N", abs(N_fit-N_true)/N_true < 1e-6, f"N_fit={N_fit:.4e}")
check("Mott-Schottky fit recovers Vbi", abs(Vbi_fit-Vbi_true) < 1e-6, f"Vbi_fit={Vbi_fit:.4f} V")

# same-code-path check: hf_cv_curve with Cox->infinity behaves like a
# bare Schottky junction (Q_sc alone sets 1/C^2 linear in V, matching Eq 13)
Cox_huge_pa = 1e6  # F/cm^2, effectively infinite oxide capacitance
phi_s_schottky = cgv.surface_potential_from_bias(np.array([VFB2+1.0]), Cox_huge_pa, VFB2, N_true, 'p')
Cs_only = cgv.semiconductor_capacitance(phi_s_schottky, N_true, 'p')
C_series_huge = cgv.total_capacitance_series_per_area(Cox_huge_pa, Cs_only)
check("series formula with Cox->inf reduces to bare Cs (same code path as MOS case)",
      abs(C_series_huge[0]-Cs_only[0])/Cs_only[0] < 1e-6,
      f"series={C_series_huge[0]:.6e}, Cs={Cs_only[0]:.6e} F/cm^2")

print("\n== 11. conductance method: single-level peak and Dit extraction ==")
Dit_true = 5e11  # cm^-2 eV^-1... (per unit area here, treated as areal density)
tau_true = 2e-5  # s
omega_grid = 2*np.pi*np.logspace(2, 7, 200)
Gpw = cgv.conductance_lorentzian(omega_grid, Dit_true, tau_true, area)
Dit_fit, tau_fit = cgv.fit_dit_from_peak(omega_grid, Gpw, area)
check("Dit recovered from noiseless single-level peak", abs(Dit_fit-Dit_true)/Dit_true < 0.02,
      f"Dit_fit={Dit_fit:.3e}, true={Dit_true:.3e}")
check("tau_it recovered from peak position", abs(tau_fit-tau_true)/tau_true < 0.02,
      f"tau_fit={tau_fit:.3e} s, true={tau_true:.3e} s")
check("peak occurs at omega*tau=1", abs(omega_grid[np.argmax(Gpw)]*tau_true - 1) < 0.05)

print("\n== 12. admittance_to_parallel: series correction identity ==")
Cm_test, Gm_test = cgv.synthetic_conductance_sweep(
    np.logspace(2, 6, 80), Dit_true, tau_true, area, tox, noise_frac=0.0)
omega_test = 2*np.pi*np.logspace(2, 6, 80)
Gpw_recovered = cgv.admittance_to_parallel(Cm_test, Gm_test, omega_test, Cox_total)
Dit_recovered, tau_recovered = cgv.fit_dit_from_peak(omega_test, Gpw_recovered, area)
check("Dit recovered through the full measure->transform->fit pipeline within 10%",
      abs(Dit_recovered-Dit_true)/Dit_true < 0.10,
      f"Dit_recovered={Dit_recovered:.3e}, true={Dit_true:.3e}")

print("\n== 13. series resistance extraction and correction ==")
Rs_true = 30.0  # ohm
f_test = 1e6
omega_1M = 2*np.pi*f_test
Cm_acc_true = Cox_total  # deep accumulation: DUT looks like Cox alone (+Rs)
# The Rs-extraction formula (Eq. 18) assumes the DUT itself is purely
# capacitive in strong accumulation (Barnes: "in good diodes G -> 0"),
# with the entire measured loss coming from Rs - simulate exactly that:
# pure Cox in series with Rs, converted back to a measured admittance.
Zmeas = Rs_true + 1/(1j*omega_1M*Cm_acc_true)
Ymeas = 1/Zmeas
Gm_meas, Cm_meas = Ymeas.real, Ymeas.imag/omega_1M
Rs_extracted = cgv.extract_series_resistance(Cm_meas, Gm_meas, omega_1M)
check("series resistance recovered exactly for a pure Cox+Rs DUT",
      abs(Rs_extracted-Rs_true)/Rs_true < 0.01, f"Rs_extracted={Rs_extracted:.3f} ohm, true={Rs_true:.1f}")
C_adj, G_adj = cgv.series_resistance_correction(Cm_meas, Gm_meas, omega_1M, Rs_extracted)
check("series-resistance-corrected capacitance recovers true Cox within 15%",
      abs(C_adj-Cm_acc_true)/Cm_acc_true < 0.15, f"C_adj={C_adj:.4e}, true Cox={Cm_acc_true:.4e} F")

print("\n== 14. mobile-charge hysteresis ==")
dVFB_true = -0.4
Nm = cgv.mobile_charge_from_hysteresis(dVFB_true, Cox_total, area)
Nm_expected = -Cox_total*dVFB_true/(cgv.Q*area)
check("mobile_charge_from_hysteresis matches definition exactly",
      abs(Nm-Nm_expected) < 1e-6, f"Nm={Nm:.3e} cm^-2")
check("mobile_charge_from_hysteresis: negative dVFB (ions migrate to Si "
      "interface on reverse sweep) gives a positive ion density",
      Nm > 0, f"Nm={Nm:.3e} cm^-2 for dVFB_true={dVFB_true:.2f} V")

print(f"\n{'='*60}\n{ok} PASSED, {fail} FAILED\n{'='*60}")
if fail:
    raise SystemExit(1)
