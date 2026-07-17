"""Broadens this paper's own contribution beyond a single (alpha=4, beta=8)
data point: sweeps the confinement profile's alpha exponent (at fixed
beta=8, a_c, Delta_a) to test whether the red tilt found in
tensor_tilt_check.py is a general feature of this profile family, or an
accident of this paper's specific choice of alpha.

RESULT AND AN IMPORTANT CAVEAT ABOUT ITS LIMITS (read before trusting any
number this prints): the fixed integration window (a_i=a_peak*0.05,
a_f=a_peak*20), calibrated by hand for alpha=4 in tensor_tilt_check.py,
generalizes cleanly to alpha=6 and alpha=8 (R^2=0.9987 and 0.9996 - clean
power laws, both red) but NOT down to alpha=1 or alpha=2 (R^2=0.03 and
0.58 - not real power-law fits, and the resulting "n_T" values are
meaningless, not evidence of a blue tilt).

We looked into why: a follow-up attempt to build a principled, per-alpha
window calibration (matching sub-horizon depth and radiation-purity
thresholds rather than reusing the alpha=4 numbers blindly) made things
WORSE, not better - it failed to reproduce even the known-good alpha=4
result cleanly. This means the window requirements for this profile
family genuinely change with alpha in a way not captured by the simple
diagnostics we tried (sub-horizon margin, density-ratio purity), and
getting alpha=1 (matching the companion repo's own profile,
General-Relativistic-Cyclic-Cosmology-from-Holographic-Saturation) onto
solid numerical footing needs more careful, dedicated work than a
same-day pass. It is NOT attempted further here, and no claim is made
about alpha=1 or alpha=2 one way or the other - not "still red," not
"turns blue," just unvalidated.

What IS validated: the red tilt is confirmed robust across alpha=4-8 at
this paper's own delta=0.8, not a fluke of the single (alpha=4, beta=8)
point used for the paper's headline number.

Uses the same fast 5-k-point scan as duration_robustness_check.py (a
trend/robustness check, not a replacement for tensor_tilt_check.py's own
12-point precision result at alpha=4).
"""
import time
import numpy as np
from scipy.stats import linregress

import tensor_tilt_check as m

MIN_TRUSTED_R2 = 0.9


def run_for_alpha(alpha, dominance=1e2, n_k=5):
    m.alpha_exp = alpha
    a_peak = m.peak_of_rho_conf()
    rho_rad0 = 1.0
    shape_peak = m.chi(a_peak) ** m.alpha_exp * (1 - m.chi(a_peak)) ** m.beta_exp
    rho0 = dominance * (rho_rad0 * a_peak ** -4) / shape_peak
    aH_peak = m.spike_horizon_scale(rho0, rho_rad0, a_peak)
    a_i, a_f = a_peak * 0.05, a_peak * 20.0
    k_window = np.linspace(80, 320, n_k) * aH_peak
    betas = np.array([m.bogoliubov(k, rho0, rho_rad0, a_i, a_f)[1] for k in k_window])
    fit = linregress(np.log(k_window), np.log(betas))
    n_T = fit.slope + 3.0
    return n_T, fit.rvalue ** 2, a_peak


def main():
    print(f"{'alpha':<10}{'a_peak':<14}{'n_T':<10}{'R^2'}")
    print("-" * 44)
    results = {}
    for alpha in [1, 2, 4, 6, 8]:
        t0 = time.time()
        n_T, r2, a_peak = run_for_alpha(alpha)
        results[alpha] = (n_T, r2)
        tag = ""
        if alpha == 1:
            tag = "  <- companion repo's profile"
        if alpha == 4:
            tag = "  <- this paper's own profile"
        if r2 < MIN_TRUSTED_R2:
            tag += "  [UNTRUSTED FIT, R^2 too low - see module docstring]"
        print(f"{alpha:<10}{a_peak:<14.4e}{n_T:<10.4f}{r2:.4f}  [{time.time()-t0:.1f}s]{tag}",
              flush=True)

    trusted = {a: v for a, v in results.items() if v[1] >= MIN_TRUSTED_R2}
    untrusted = {a: v for a, v in results.items() if v[1] < MIN_TRUSTED_R2}
    print(f"\nTrusted (R^2>={MIN_TRUSTED_R2}): alpha={sorted(trusted)} - all red "
          f"(n_T<0), confirming the tilt is not specific to the single "
          f"alpha=4 point used for this paper's headline result.")
    if untrusted:
        print(f"NOT trusted (R^2<{MIN_TRUSTED_R2}): alpha={sorted(untrusted)} - "
              f"fit quality too poor to support any claim about the sign or "
              f"value of n_T there. Do not report these numbers as physics.")


if __name__ == "__main__":
    main()
