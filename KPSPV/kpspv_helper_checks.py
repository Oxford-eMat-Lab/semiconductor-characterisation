"""Analytic-limit checks for KPSPV/kpspv_helper.py.

Not part of the repo layout (the other techniques have no test file) - run
it from inside KPSPV/ with:  python3 kpspv_helper_checks.py
"""
import numpy as np, kpspv_helper as kp

ok = fail = 0
def check(name, cond, detail=""):
    global ok, fail
    if cond: ok += 1; print(f"  PASS  {name}   {detail}")
    else:    fail += 1; print(f"  FAIL  {name}   {detail}")

print("== 1. CPD and calibration ==")
check("CPD = 0 for equal work functions", kp.cpd_from_work_functions(5.1, 5.1) == 0.0)
c = kp.cpd_from_work_functions(4.28, 5.10)
check("Au probe on Al sample gives -0.82 V", abs(c + 0.82) < 1e-9, f"CPD={c:.4f} V")
phi_probe = kp.calibrate_probe(kp.cpd_from_work_functions(4.60, 5.10), 4.60)
check("calibrate_probe round-trips", abs(phi_probe - 5.10) < 1e-12, f"Phi_probe={phi_probe:.6f} eV")
back = kp.work_function_from_cpd(c, 5.10)
check("work_function_from_cpd inverts", abs(back - 4.28) < 1e-12, f"Phi_s={back:.6f} eV")

print("\n== 2. the off-null fit ==")
Vb = np.linspace(-1.5, 0.5, 21)
cpd_true = -0.412
y = kp.offnull_amplitude(Vb, cpd_true, 3.0)
r = kp.fit_null_point(Vb, y)
check("noiseless sweep recovers CPD exactly", abs(r.cpd_V - cpd_true) < 1e-10, f"{r.cpd_V:.9f} V")
check("gradient recovered", abs(r.gradient - 3.0) < 1e-10, f"k={r.gradient:.6f}")
y2 = kp.synthetic_null_sweep(Vb, cpd_true, 3.0, noise_mV=2.0, spacing_drift=0.0, seed=0)
r2 = kp.fit_null_point(Vb, y2)
check("noisy sweep within 3 sigma", abs(r2.cpd_V - cpd_true) < 3*r2.cpd_err, f"{r2.cpd_V*1e3:.2f}+/-{r2.cpd_err*1e3:.2f} mV")
y3 = kp.synthetic_null_sweep(Vb, cpd_true, 3.0, noise_mV=0.0, spacing_drift=0.25, seed=0)
r3 = kp.fit_null_point(Vb, y3)
check("a steady but different spacing changes k, not the crossing",
      abs(kp.fit_null_point(Vb, kp.offnull_amplitude(Vb, cpd_true, 8.0)).cpd_V - cpd_true) < 1e-10)
check("drift DURING a sweep biases the crossing", abs(r3.cpd_V - cpd_true) > 0.01, f"biased to {r3.cpd_V*1e3:.2f} mV from {cpd_true*1e3:.1f} mV, R^2={r3.r_squared:.6f}")

print("\n== 3. the vibrating capacitor ==")
t = np.linspace(0, 1/80.0, 4001)
i_null = kp.kelvin_current(t, -cpd_true, cpd_true)
check("current vanishes identically at null", np.allclose(i_null, 0.0), f"max|i|={np.abs(i_null).max():.2e} A")
i_off = kp.kelvin_current(t, -cpd_true + 0.1, cpd_true)
check("off-null current is nonzero and AC", np.abs(i_off).max() > 0 and abs(np.mean(i_off)) < 1e-3*np.abs(i_off).max())
C = kp.kelvin_capacitance(t, d1_um=0.0)
check("static capacitance = eps0 A/d", np.allclose(C, kp.EPS0*3.14e-2/(200e-4)))
check("modulation index", kp.modulation_index(200., 40.) == 0.2)

print("\n== 4. bulk semiconductor ==")
Vt = kp.thermal_voltage()
check("thermal voltage ~25.85 mV", abs(Vt-0.025852) < 1e-5, f"{Vt*1e3:.4f} mV")
n_b, p_b = kp.bulk_carrier_densities(5e15, 'n')
check("np = ni^2 in equilibrium", abs(n_b*p_b/kp.NI_SI**2 - 1) < 1e-9, f"n_b={n_b:.3e}, p_b={p_b:.3e}")
check("n_b = N_D for N_D >> ni", abs(n_b/5e15 - 1) < 1e-9)
phiF = kp.fermi_potential(5e15, 'n')
check("phi_F positive n-type, ~0.34 V", 0.33 < phiF < 0.35, f"{phiF:.4f} V")
check("phi_F sign flips for p-type", kp.fermi_potential(5e15,'p') == -phiF)
Ws = kp.work_function_semiconductor(5e15, 'n')
check("n-Si 5e15 work function 4.2-4.4 eV", 4.2 < Ws < 4.4, f"Phi_s={Ws:.4f} eV (chi=4.05; chi=4.10 would give {Ws+0.05:.4f})")

print("\n== 5. space charge: analytic limits ==")
check("Q_sc = 0 at flat band", kp.space_charge_density(0.0, 5e15) == 0.0)
for phi in (-0.2, -0.35, -0.5):
    Q = kp.space_charge_density(phi, 5e15, 'n')
    Qd = np.sqrt(2*kp.K_SI*kp.EPS0*kp.Q*5e15*(abs(phi)-Vt))/kp.Q   # exact limit carries -kT/q
    Qnaive = np.sqrt(2*kp.K_SI*kp.EPS0*kp.Q*5e15*abs(phi))/kp.Q
    check(f"depletion limit (phi_s - kT/q) at phi_s={phi} V", abs(Q/Qd - 1) < 5e-3,
          f"full={Q:.4e}, depl(phi-kT/q)={Qd:.4e}, ratio={Q/Qd:.5f}; naive depl overestimates by {100*(Qnaive/Q-1):.1f}%")
check("Q_sc > 0 for phi_s < 0 (n-type depletion)", kp.space_charge_density(-0.3, 5e15) > 0)
check("Q_sc < 0 for phi_s > 0 (n-type accumulation)", kp.space_charge_density(+0.3, 5e15) < 0)
# Schroder's closed form 9.07e-7 N^2/(Ks N_dop)
N = 1e11; NA = 1e16
phi_schroder = kp.band_bending_depletion(N, NA, K_s=11.7)
check("Schroder worked example: 1e11 cm^-2 on 1e16 -> 0.077 V", abs(phi_schroder-0.077) < 0.001, f"{phi_schroder:.4f} V")
check("prefactor is 9.05e-7", abs(kp.Q/(2*kp.EPS0)/1e-7 - 9.05) < 0.05, f"{kp.Q/(2*kp.EPS0):.4e}")
# round trip: depletion phi -> Q -> phi
Qtest = kp.space_charge_density(-0.30, 1e16, 'n')
check("Q_sc -> band_bending_depletion recovers phi_s - kT/q",
      abs(kp.band_bending_depletion(Qtest, 1e16) - (0.30 - Vt)) < 1e-3,
      f"{kp.band_bending_depletion(Qtest,1e16):.5f} V vs phi_s-kT/q = {0.30-Vt:.5f} V")
W = kp.depletion_width_um(0.3, 5e15)
check("depletion width ~0.28 um at 0.3 V, 5e15", 0.2 < W < 0.4, f"W={W:.4f} um")
check("W -> 0 as phi_s -> 0", kp.depletion_width_um(0.0, 5e15) == 0.0)

print("\n== 6. dielectric term ==")
check("V_i = 0 for zero charge", kp.insulator_potential(0.0, 100.0) == 0.0)
check("V_i = 0 for centroid at interface", kp.insulator_potential(1e12, 0.0) == 0.0)
Vi = kp.insulator_potential(1e12, 100.0, kp.K_SIO2)
check("1e12 cm^-2 at 100 nm in SiO2 -> ~4.6 V", 4.0 < Vi < 5.2, f"{Vi:.4f} V")
check("V_i linear in Q_f", abs(kp.insulator_potential(2e12,100.)/Vi - 2) < 1e-12)
check("V_i sign follows charge sign", kp.insulator_potential(-1e12,100.) < 0)

print("\n== 7. interface states ==")
q0 = kp.interface_charge(0.0, 5e15, 'n', dit_midgap=0.0, dit_edge=0.0)
check("Q_it = 0 for zero Dit", q0 == 0.0)
qa = kp.interface_charge(+0.3, 5e15, 'n')   # accumulation: E_F near E_C -> acceptors filled
qd = kp.interface_charge(-0.6, 5e15, 'n')   # inversion: E_F near E_V -> donors emptied
check("Q_it more negative in accumulation than inversion", qa < qd, f"acc={qa:.3e}, inv={qd:.3e}")
check("Q_it changes sign across the gap", qa < 0 < qd, f"acc={qa:.3e}, inv={qd:.3e}")

print("\n== 8. charge balance ==")
phi_ideal = kp.solve_surface_potential(0.0, 5e15, 'n', dit_midgap=1e-30, dit_edge=1e-30)
check("phi_s = 0 for zero charge AND zero Dit", abs(phi_ideal) < 1e-6, f"phi_s={phi_ideal*1e6:.3f} uV")
phi0 = kp.solve_surface_potential(0.0, 5e15, 'n')
qit0 = kp.interface_charge(phi0, 5e15, 'n'); qsc0 = kp.space_charge_density(phi0, 5e15, 'n')
check("with Dit, an UNCHARGED dielectric still bends the bands", -0.15 < phi0 < -0.005,
      f"phi_s={phi0*1e3:.2f} mV; Q_it={qit0:.3e} balanced by Q_sc={qsc0:.3e} cm^-2")
phi_pos = kp.solve_surface_potential(+2e11, 5e15, 'n')
phi_neg = kp.solve_surface_potential(-2e11, 5e15, 'n')
check("positive Qf -> accumulation (phi_s>0) on n-type", phi_pos > 0, f"{phi_pos:.4f} V")
check("negative Qf -> depletion/inversion (phi_s<0) on n-type", phi_neg < 0, f"{phi_neg:.4f} V")
check("regime names", (kp.surface_regime(phi_pos)=='accumulation'
                       and kp.surface_regime(phi_neg) in ('depletion','inversion')
                       and kp.surface_regime(0.0)=='flat band'),
      f"{kp.surface_regime(phi_pos)} / {kp.surface_regime(phi_neg)}")
res = 2e11 + kp.space_charge_density(phi_pos,5e15) + kp.interface_charge(phi_pos,5e15)
check("neutrality residual is zero at the solution", abs(res) < 1e4, f"residual={res:.3e} cm^-2")
phi_mono = [kp.solve_surface_potential(q, 5e15, 'n') for q in np.linspace(-4e11, 4e11, 9)]
check("phi_s increases monotonically with Qf", np.all(np.diff(phi_mono) > 0))

print("\n== 9. light and SPV ==")
cd_, cl_, spv_, pd_, pl_ = kp.cpd_dark_and_light(-2e11, 1e15, 5e15, 'n', xc_nm=5.0)
check("SPV = CPD_dark - CPD_light", abs(spv_ - (cd_-cl_)) < 1e-15)
check("SPV = -(phi_dark - phi_light)", abs(spv_ + (pd_-pl_)) < 1e-12, f"SPV={spv_*1e3:.2f} mV")
check("illumination flattens the bands", abs(pl_) < abs(pd_), f"dark={pd_:.4f} V, light={pl_:.4f} V")
check("negative Qf on n-type gives POSITIVE SPV", spv_ > 0, f"SPV={spv_*1e3:.2f} mV")
_,_,spv_acc,_,_ = kp.cpd_dark_and_light(+2e11, 1e15, 5e15, 'n', xc_nm=5.0)
check("positive Qf on n-type gives NEGATIVE SPV", spv_acc < 0, f"SPV={spv_acc*1e3:.2f} mV")
_,_,spv0,_,_ = kp.cpd_dark_and_light(-2e11, 0.0, 5e15, 'n')
check("SPV = 0 with no injection", abs(spv0) < 1e-12, f"{spv0:.2e} V")
dns = np.logspace(11, 18, 8)
spvs = np.array([kp.cpd_dark_and_light(-2e11, d, 5e15,'n')[2] for d in dns])
check("SPV grows monotonically with injection", np.all(np.diff(spvs) > -1e-9))
check("SPV saturates at -phi_s(dark)", abs(spvs[-1] + pd_) < 0.02, f"sat={spvs[-1]:.4f} V, -phi_dark={-pd_:.4f} V")
_,_,_,pdk,plt_ = kp.cpd_dark_and_light(-2e11, 1e18, 5e15,'n')
check("high injection drives flat band", abs(plt_) < 0.03, f"phi_light={plt_*1e3:.2f} mV")

print("\n== 10. lifetime and S ==")
check("tau_eff -> tau_bulk as S -> 0", abs(kp.effective_lifetime_us(1000., 0.0, 200.)-1000.) < 1e-6)
t_hi = kp.effective_lifetime_us(1e9, 1e7, 200.)
check("tau_eff -> W/2S as S -> inf", abs(t_hi - 1e6*200e-4/(2*1e7)) < 1e-6*t_hi, f"{t_hi:.4e} us")
check("tau_eff decreases with S", kp.effective_lifetime_us(1000.,10.,200.) > kp.effective_lifetime_us(1000.,1000.,200.))
S_flat = kp.surface_recombination_velocity(0.0, 5e15,'n', delta_n_cm3=1e14)
S_acc  = kp.surface_recombination_velocity(0.4, 5e15,'n', delta_n_cm3=1e14)
check("S falls when bands are bent (field effect)", S_acc < S_flat, f"flat={S_flat:.3g}, bent={S_acc:.3g} cm/s")
check("S -> 0 as Dit -> 0", kp.surface_recombination_velocity(0.0,5e15,'n',dit_midgap=1e-30,dit_edge=1e-30) < 1e-6)

print("\n== 11. Goodman diffusion length ==")
wl = np.array([800,850,900,950,1000,1030,1060])
a = kp.alpha_silicon(wl)
check("alpha decreasing with wavelength", np.all(np.diff(a) < 0))
Vp = kp.spv_constant_flux(a, 250.0, C2=1.0)
L, Lerr, _ = kp.fit_diffusion_length(a, Vp)
check("noiseless Goodman fit recovers L_n exactly", abs(L-250.) < 1e-6, f"L_n={L:.6f} um")
Vpn = kp.synthetic_spv_wavelength_scan(wl, 250.0, noise_frac=0.02, seed=3)
Ln, Lne, _ = kp.fit_diffusion_length(a, Vpn)
check("noisy Goodman fit within 3 sigma", abs(Ln-250.) < 3*Lne, f"L_n={Ln:.1f}+/-{Lne:.1f} um")
check("uncertainty is large (intercept extrapolation)", Lne/Ln > 0.02, f"{100*Lne/Ln:.1f}% error")

print("\n== 12. charge fluctuations ==")
f = lambda q: kp.solve_surface_potential(q, 5e15, 'n')
check("sigma=0 reproduces the point value", kp.charge_fluctuation_average(1e11,0.0,f) == f(1e11))
smeared = [kp.charge_fluctuation_average(q, 2e11, f) for q in np.linspace(-4e11,4e11,9)]
sharp   = [f(q) for q in np.linspace(-4e11,4e11,9)]
check("fluctuations reduce the phi_s swing", (max(smeared)-min(smeared)) < (max(sharp)-min(sharp)),
      f"smeared={max(smeared)-min(smeared):.3f} V, sharp={max(sharp)-min(sharp):.3f} V")

print(f"\n{'='*60}\n{ok} passed, {fail} failed\n{'='*60}")
raise SystemExit(1 if fail else 0)
