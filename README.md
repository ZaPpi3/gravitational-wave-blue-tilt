# gravitational-wave-blue-tilt

Analytical derivations and tensor perturbation tracking pipelines demonstrating the organic generation of a highly blue-tilted primordial gravitational wave spectrum from a transient hyper-stiff brane rupture phase.

## Overview

In rupture-driven cosmology, the topological decoupling of an emergent brane universe triggers a transient confinement-release energy spike where the effective equation of state peaks near $w_{\text{conf}} \approx 80$ prior to recombination. 

Because the background energy density dilutes exponentially faster than standard radiation ($\rho_{\text{conf}} \propto a^{-243}$), sub-horizon tensor perturbations experience massive relative amplification. This mechanism natively generates a highly blue-tilted primordial tensor spectrum ($n_T \simeq +2.0$) without violating the Null Energy Condition or invoking fine-tuned scalar fields.

The resulting power-law spectrum naturally accounts for the enhanced stochastic gravitational wave background reported by recent pulsar timing arrays, while leaving a distinct signature uniquely testable by future space-based interferometers like LISA and DECIGO.

## Repository Contents

* `main.tex`: Full publication-ready LaTeX source file.
* `main.pdf`: Compiled, high-resolution rapid-communication Letter layout.
* `code/make_figure_gw.py`: Optimized Python implementation script utilizing Matplotlib to calculate and generate the tensor perturbation tracking curves.
* `figures/`: Output directory hosting the vector plot showing the blue-tilted power-law scaling and target sensor sensitivity windows.

## Reference Parameters

The tracking architecture evaluates the fiducial cosmology using the following parameters:
$$\alpha = 4, \quad \beta = 8, \quad w_{\max} = 80.0, \quad f_{\text{anchor}} = 10^{-9}\,\text{Hz}, \quad \Omega_{\text{PTA}} = 10^{-14}$$

---
**ORCID:** [0009-0009-8933-857X](https://orcid.org)
