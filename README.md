# Gravitational-Wave Blue Tilt

Analytical derivations and a numerical Mukhanov-Sasaki integration for the tensor spectrum sourced by a transient hyper-stiff brane rupture phase.

**Correction (2026-07-17):** this repository originally claimed a blue-tilted tensor spectrum ($n_T \simeq +2.0$) as a possible explanation for PTA observations. That number came from a closed-form formula valid only for an eternal, unchanging equation of state - not the finite-duration transient this paper actually describes - and was presented alongside a Python code listing that looked like a real integration but never actually called a solver. An honest Mukhanov-Sasaki integration of this paper's own confinement profile (`code/tensor_tilt_check.py`) gives the opposite result: a **red** tilt, $n_T \approx -0.9$. No PTA-band signal is predicted by this mechanism. This is consistent with an independent finding on a related profile in the companion repository `General-Relativistic-Cyclic-Cosmology-from-Holographic-Saturation` ($n_T\approx-1$). See `main.tex` for the full corrected derivation.

## Overview
This repository studies a model where a transient hyper-stiff phase ($w_{\text{conf}} \approx 80$ in the paper's own $w(a)=w_{\text{conf}}\chi(a)$ ansatz, evaluated far from the density peak - see correction below) sources tensor perturbations during a confinement-release energy spike in rupture-driven brane cosmology. The background solution and its rapid dilution rate are unaffected by the correction above; the tensor-spectrum prediction and PTA/LISA/DECIGO observational claims are withdrawn and replaced by the computed red tilt.

**Second correction (same date):** an earlier version of this fix also claimed the *self-consistent* effective equation of state (derived from the continuity equation, not the ansatz above) "peaks near w~80 close to the density peak." That was wrong too - it was a boundary artifact of an arbitrarily narrow numerical search window. The self-consistent w(a) is exactly $-1$ *at* the density peak (a general fact about local maxima of any smooth density profile) and only grows large far into the tail where the confinement density itself has become negligible. This does not affect the computed tensor tilt, which integrates the full w(a) trajectory directly rather than relying on any single "peak w" figure - but the descriptive text was wrong and has been corrected in `main.tex` and `code/tensor_tilt_check.py`.

**Robustness addendum (same date):** the red tilt is confirmed stable against the confinement profile's shape exponent alpha (tested at alpha=4,6,8, all clean power-law fits) and against widening the transition duration by up to 100x (still converges to n_T~-1). An attempt to extend the alpha-robustness check down to alpha=1 - the exponent used by the companion repo `General-Relativistic-Cyclic-Cosmology-from-Holographic-Saturation` - did not produce reliable fits with either the original or a more carefully recalibrated integration window; we report this as an open question rather than force a number either way. See `main.tex` for the full discussion.

## Repository Contents
- `main.tex`: LaTeX source.
- `main.pdf`: Compiled manuscript.
- `code/tensor_tilt_check.py`: The real Mukhanov-Sasaki integration - builds the confinement + radiation background, derives the confinement fluid's effective equation of state from the continuity equation (not assumed), integrates the tensor mode equation with adaptive Runge-Kutta, extracts Bogoliubov coefficients, and regenerates `figures/figure_gw.png`.
- `code/duration_robustness_check.py`: Sanity check - does widening the confinement transition (Delta_a) push the tilt back toward blue, as the eternal-w formula would suggest? No: widening by up to 100x still converges to n_T~-1 (R^2->1.0000), showing the red tilt is a robust attractor of this profile family, not a fragile artifact of the paper's specific Delta_a.
- `code/profile_family_sweep.py`: Sanity check across the profile's shape exponent alpha (beta=8 fixed). Red tilt confirmed robust for alpha=4,6,8 (R^2=0.955-0.9996). An attempt to extend this down to alpha=1-2 - which would bridge directly to the companion repo's own profile - did **not** produce reliable fits (R^2=0.03-0.58) despite a follow-up calibration attempt; this is reported honestly as unresolved, not swept under the rug.
- `figures/figure_gw.png`: Computed $|\beta_k|^2$ spectrum and power-law fit.

## Reference Simulation Parameters
$\alpha = 4$, $\beta = 8$, $a_c = 5\times10^{-4}$, $\Delta a = 4\times10^{-4}$. ($w_{\max}\approx80$ refers to the paper's separate $w(a)=w_{\text{conf}}\chi(a)$ ansatz, not to the self-consistent w(a) at the density peak, which is exactly $-1$ there - see `main.tex` Sec. II.)

## Installation & Reproducibility
```bash
git clone https://github.com/ZaPpi3/gravitational-wave-blue-tilt
cd gravitational-wave-blue-tilt
pip install numpy scipy matplotlib
python code/tensor_tilt_check.py
```
This is slow: each $k$ takes roughly 15-70+ seconds to integrate (the solver must resolve both the confinement spike and many oscillation periods at high $k$), so the default 12-point, two-amplitude scan takes on the order of 15-20 minutes.

## Citation
Please cite the accompanying manuscript.
```bibtex
@article{jarvis2026tensorpert,
  title={Tensor Perturbations from a Transient Hyper-Stiff Brane Rupture Phase},
  author={Jarvis, Paul},
  year={2026}
}
```

## License
MIT License.
