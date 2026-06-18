# Gravitational-Wave Blue Tilt

Analytical derivations and numerical tracking pipelines for a highly blue-tilted primordial gravitational wave spectrum from a transient hyper-stiff brane rupture phase.

## Overview
This repository provides the code and manuscript associated with a model where a transient hyper-stiff phase (peak $w_{\text{conf}} \approx 80$) generates a blue-tilted tensor spectrum ($n_T \simeq +2.0$). This mechanism offers a potential explanation for pulsar timing array observations.

## Repository Contents
- `manuscript/main.tex`: LaTeX source.
- `code/make_figure_gw.py`: Python script for numerical simulation and plotting.
- `figures/figure_gw.png`: Output visualization.

## Reference Simulation Parameters
The code uses the following parameters for the fiducial model:
$\alpha = 4$, $\beta = 8$, $w_{\max} = 80.0$.

## Installation & Reproducibility
Clone the repo and run the numerical simulation:
```bash
git clone https://github.com
cd gravitational-wave-blue-tilt
pip install numpy scipy matplotlib
python code/make_figure_gw.py
```

## Citation
Please cite the accompanying manuscript.
```bibtex
@article{jarvis2026bluetilted,
  title={Blue-Tilted Primordial Gravitational Waves from a Transient Hyper-Stiff Brane Rupture Phase},
  author={Jarvis, Paul},
  year={2026}
}
```

## License
MIT License.
