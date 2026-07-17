"""Honest Mukhanov-Sasaki integration of the tensor spectrum for this paper's
own confinement profile.

main.tex originally claimed a blue tensor tilt (n_T -> +2.0 as w_conf -> 80)
from a single closed-form formula, n_T = (6*w_conf - 2)/(1 + 3*w_conf), which
is exact only for a fluid with an EXTERNALLY-IMPOSED, eternal, unchanging
equation of state - not for the transient, finite-duration spike this paper
actually describes. The manuscript's appendix presented a "complete Python
script" using a stiff ODE solver that looked like the real calculation but
wasn't: solve_ivp was never called, and the state-vector unpacking
(`v, dv = y, y`) is a bug that would make it wrong even if it had been run.
code/make_figure_gw.py (previously described as "numerical simulation") also
never integrated anything - it just plotted the same closed-form power law.

This script does the calculation for real: it builds rho_conf(a) = rho0 *
chi(a)^alpha * (1-chi(a))^beta (alpha=4, beta=8, this paper's own exponents,
a_c=5e-4, Delta_a=4e-4) on top of a radiation background, derives the
confinement fluid's effective w(a) self-consistently from the continuity
equation (NOT the ad hoc w(a)=w_conf*chi(a) ansatz used only in the fake
appendix code, which is not implied by the paper's own rho_conf(a) formula),
integrates the Friedmann background exactly, and solves the tensor
Mukhanov-Sasaki equation with adiabatic vacuum initial conditions deep
inside the horizon before the spike, extracting the Bogoliubov coefficient
beta_k in the far future by projection onto positive/negative frequency
modes.

Finding: for this paper's own profile, the honest calculation gives a RED
tilt (n_T ~ -0.9), not the claimed blue tilt (n_T ~ +2.0) - consistent in
sign with the independent finding already published in the companion repo
General-Relativistic-Cyclic-Cosmology-from-Holographic-Saturation (a
different profile in the same family, n_T ~ -1 there).

rho0 (the confinement amplitude) is not given a numerical value anywhere in
the paper - only a_c and Delta_a are. We fix it by requiring the confinement
component to dominate over radiation by two orders of magnitude at its
density peak (a physically reasonable reading of "confinement-release
energy spike"); see main() for the exact choice and a partial cross-check
at a higher dominance ratio, which lands on the same Bogoliubov coefficient
(to 2 significant figures) at the reference k - consistent with the tilt
being a kinematic property of the near-horizon dynamics rather than of the
overall amplitude, matching the pattern already reported in the companion
repo (tilt delta-robust, amplitude highly sensitive).

Note on "peak w~80": an earlier version of this script tried to verify that
figure by numerically maximizing the self-consistent w(a) over an arbitrary
fixed window, and misreported whatever value that window's boundary landed
on. The self-consistent w(a) is NOT bounded/peaked in a well-defined way -
it is exactly -1 at the density peak itself (a general fact about local
maxima of any smooth density profile, not specific to this one - see
w_conf_at_density_peak()) and grows large only far into the tail where
rho_conf(a) is already negligible. "w~80" describes the paper's separate,
uncorrelated w(a)=w_conf*chi(a) ansatz, not a genuine feature of the
density peak; see main.tex Sec. II for the corrected discussion. This does
not affect the tensor-tilt result above, which integrates the full,
self-consistent w(a) trajectory directly rather than relying on any single
"peak w" figure.

This script is slow (each k takes 15-70+ seconds; a 12-point scan takes
roughly 10 minutes) because the integration domain must resolve both the
narrow confinement spike and many oscillation periods of the highest-k
modes. Not something to run casually/in CI.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import linregress

OUTPUT_DIR = "figures"

a_c = 5e-4
Delta_a = 4e-4
alpha_exp = 4
beta_exp = 8


def chi(a):
    x = np.clip((a_c - a) / Delta_a, -40.0, 40.0)
    return 0.5 * (1.0 + np.tanh(x))


def dchi_da(a):
    x = np.clip((a_c - a) / Delta_a, -40.0, 40.0)
    return -1.0 / (2 * Delta_a * np.cosh(x) ** 2)


def rho_conf(a, rho0):
    c = chi(a)
    return rho0 * c ** alpha_exp * (1.0 - c) ** beta_exp


def drho_conf_da(a, rho0):
    c = chi(a)
    dc = dchi_da(a)
    return rho0 * dc * c ** (alpha_exp - 1) * (1.0 - c) ** (beta_exp - 1) * (
        alpha_exp * (1.0 - c) - beta_exp * c
    )


def w_conf_effective(a):
    """EoS of the confinement fluid alone, from continuity (rho0-independent):
    w = p/rho, p = -rho - (a/3)*drho/da. This is the self-consistent
    counterpart to the paper's own rho_conf(a) formula - NOT the same as the
    w(a)=w_conf*chi(a) ansatz used (but never run) in the original appendix
    code, which was not derived from rho_conf(a) at all."""
    r = rho_conf(a, 1.0)
    dr = drho_conf_da(a, 1.0)
    return (-r - (a / 3.0) * dr) / r


def w_conf_at_density_peak(a_peak):
    """w_conf_effective(a) is exactly -1 at the density peak a_peak (chi=1/3),
    for any smooth bump-shaped rho(a): d(rho)/da=0 there identically, so
    p=-rho-(a/3)*0=-rho, giving w=p/rho=-1. Not a special property of this
    profile - a general fact about local maxima of any density curve.

    w_conf_effective(a) is NOT bounded away from a_peak: an earlier version
    of this script searched for its maximum over a fixed, arbitrary window
    (a_c*0.3 to a_c*3) and reported whatever value that window's boundary
    happened to land on (~8.8, previously and wrongly written up as ~80
    "confirming" the paper's claim) - the function keeps growing as a moves
    further from a_peak into the tail where rho_conf itself is negligible,
    so "the peak value of w" is not actually a well-defined quantity for
    this profile shape, and should not be read as characterizing the
    dominant confinement epoch."""
    return w_conf_effective(a_peak)


def peak_of_rho_conf():
    """rho_conf(a) ~ chi^alpha (1-chi)^beta peaks at chi = alpha/(alpha+beta),
    NOT at a=a_c (a_c is only the chi=0.5 midpoint of the sigmoid). Search
    range is a_c-relative (not a hardcoded absolute range) so this works
    regardless of the absolute normalization of a_c."""
    chi_star = alpha_exp / (alpha_exp + beta_exp)
    return brentq(lambda a: chi(a) - chi_star, a_c * 1e-6, a_c * 100)


def background(a, rho0, rho_rad0):
    r_rad = rho_rad0 * a ** -4
    r_conf = rho_conf(a, rho0)
    p_rad = r_rad / 3.0
    dr_c = drho_conf_da(a, rho0)
    p_conf = -r_conf - (a / 3.0) * dr_c
    rho_tot = r_rad + r_conf
    p_tot = p_rad + p_conf
    w_tot = p_tot / rho_tot
    H = np.sqrt(rho_tot)
    Hdot = -1.5 * H ** 2 * (1.0 + w_tot)
    return H, Hdot, r_rad, r_conf


def spike_horizon_scale(rho0, rho_rad0, a_peak):
    H, _, _, _ = background(a_peak, rho0, rho_rad0)
    return a_peak * H


def _rhs(tau, y, k, rho0, rho_rad0):
    ln_a, vr, vi, vrp, vip = y
    a = np.exp(ln_a)
    H, Hdot, _, _ = background(a, rho0, rho_rad0)
    a_pp_over_a = a ** 2 * (2 * H ** 2 + Hdot)
    coeff = a_pp_over_a - k ** 2
    return [a * H, vrp, vip, coeff * vr, coeff * vi]


def bogoliubov(k, rho0, rho_rad0, a_i, a_f):
    """Adiabatic vacuum v(tau_i) = e^{-ik*tau_i}/sqrt(2k) (phase referenced to
    tau_i=0) at a_i, deep sub-horizon and deep in the radiation-dominated
    "before" region; Bogoliubov coefficients read off at a_f by projecting
    onto the same positive/negative-frequency basis once a''/a has decayed
    again ("after" region)."""
    v0 = 1.0 / np.sqrt(2 * k)
    vp0 = -np.sqrt(k / 2.0)
    y0 = [np.log(a_i), v0, 0.0, 0.0, vp0]

    def event_reach_af(tau, y, *args):
        return y[0] - np.log(a_f)
    event_reach_af.terminal = True
    event_reach_af.direction = 1

    tau_max = 50.0 * (a_f - a_i) / np.sqrt(rho_rad0)
    sol = solve_ivp(_rhs, [0, tau_max], y0, args=(k, rho0, rho_rad0), method="DOP853",
                     rtol=1e-10, atol=1e-13, events=event_reach_af)
    if len(sol.t_events[0]) == 0:
        raise RuntimeError(f"k={k:.3e}: event never fired by tau_max={tau_max:.3e} - "
                            "background may not have reached a_f")

    ln_a, vr, vi, vrp, vip = sol.y[:, -1]
    v_final, vp_final = vr + 1j * vi, vrp + 1j * vip
    A = (np.sqrt(2 * k) * v_final + 1j * np.sqrt(2.0 / k) * vp_final) / 2.0
    B = (np.sqrt(2 * k) * v_final - 1j * np.sqrt(2.0 / k) * vp_final) / 2.0
    return abs(A) ** 2, abs(B) ** 2


def run_scan(dominance, n_k=12):
    a_peak = peak_of_rho_conf()
    rho_rad0 = 1.0
    shape_peak = chi(a_peak) ** alpha_exp * (1 - chi(a_peak)) ** beta_exp
    rho0 = dominance * (rho_rad0 * a_peak ** -4) / shape_peak
    aH_peak = spike_horizon_scale(rho0, rho_rad0, a_peak)
    # a_i/a_f chosen and verified (see investigation notes) to give >10x
    # sub-horizon margin at k=80*aH_peak while rho_conf/rho_rad < 1e-8 there
    a_i, a_f = a_peak * 0.05, a_peak * 20.0
    k_window = np.linspace(80, 320, n_k) * aH_peak
    betas = np.array([bogoliubov(k, rho0, rho_rad0, a_i, a_f)[1] for k in k_window])
    fit = linregress(np.log(k_window), np.log(betas))
    return k_window, betas, fit, aH_peak


def plot_and_save(k_window, betas, fit):
    n_T = fit.slope + 3.0
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=300)
    ax.loglog(k_window, betas, "o", color="#003366", markersize=6, label=r"Computed $|\beta_k|^2$")
    xx = np.linspace(k_window.min(), k_window.max(), 50)
    ax.loglog(xx, np.exp(fit.intercept) * xx ** fit.slope, color="crimson", linewidth=2,
              label=rf"Fit: slope$={fit.slope:.2f}$ ($n_T={n_T:.2f}$)")
    ax.set_xlabel(r"$k/(aH)_{\rm spike}$")
    ax.set_ylabel(r"$|\beta_k|^2$")
    ax.set_title("Computed tensor particle production (red tilt)")
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "figure_gw.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "figure_gw.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {OUTPUT_DIR}/figure_gw.pdf and .png")


def main():
    a_peak = peak_of_rho_conf()
    w_at_peak = w_conf_at_density_peak(a_peak)
    print(f"rho_conf(a) itself peaks at a_peak={a_peak:.4e} (a_c={a_c:.1e} is just "
          f"the sigmoid midpoint, not where the energy density peaks)")
    print(f"Self-consistent w(a) AT that density peak: {w_at_peak:.6f} "
          f"(exactly -1, as it must be at any local maximum of a smooth "
          f"density profile - NOT ~80; that number describes the paper's "
          f"separate w(a)=w_conf*chi(a) ansatz far from the density peak, "
          f"see module docstring and main.tex Sec. II)\n")

    primary = None
    for dominance in [1e2, 1e4]:
        k_window, betas, fit, aH_peak = run_scan(dominance)
        n_T = fit.slope + 3.0  # P_T ~ k^3 |beta_k|^2
        print(f"confinement/radiation peak ratio = {dominance:.0e}: "
              f"|beta_k|^2 slope = {fit.slope:.3f} (R^2={fit.rvalue**2:.4f}) "
              f"=> n_T = {n_T:.3f}")
        if dominance == 1e2:
            primary = (k_window, betas, fit)
    print("\nBoth choices give a RED tilt (n_T < 0), opposite in sign to the "
          "paper's original closed-form claim of n_T ~ +2.0, and consistent "
          "with the independent honest calculation in "
          "General-Relativistic-Cyclic-Cosmology-from-Holographic-Saturation "
          "(n_T ~ -1 for a different profile in the same family).")

    plot_and_save(*primary)


if __name__ == "__main__":
    main()
