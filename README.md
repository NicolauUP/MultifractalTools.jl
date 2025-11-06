# MultifractalTools

[![Build Status](https://github.com/NicolauUP/MultifractalTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/NicolauUP/MultifractalTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package for performing multifractal analysis on 1D and 2D data.

## 🚀 Features

* Calculates multifractal spectra: $\tau(q)$, $\alpha(q)$, and $f(\alpha)$.
* Supports 1D vectors and 2D matrices.
* Handles zero-measure bins for robust numerical stability.
* Interactive plotting (`plot_to_fit`) via `GLMakie` to visually find scaling regions.
* Automated plotting for the results;

## 📦 Installation

This package is currently being registered. Once available, you can install it via the Julia package manager:

```julia
using Pkg
Pkg.add("MultifractalTools")
```
In the meantime, you can install it directly from GitHub:
```julia
using Pkg
Pkg.add(url="[https://github.com/NicolauUP/MultifractalTools.jl](https://github.com/NicolauUP/MultifractalTools.jl)")
```

Quick Start

Here is a minimal example of the main workflow:
```julia

using MultifractalTools

# 1. Create sample data (1D or 2D)
data = rand(128, 128)

# 2. Get q values
qs = obtain_qs(-10, 10, 21)

# 3. Compute scaling quantities
scaling_data = compute_scaling_quantities(data, qs)

# 4. Compute the spectrum (e.g., fitting from l-index 2 to 6)
spec = compute_spectrum(scaling_data, qs, 2, 6)

# 5. Plot the results
# (Requires 'using GLMakie')
plot_spectrum(spec, which=:Spectrum)
```
For an interactive plot to help you choose the fitting range (step 4), see plot_to_fit.

📖 Documentation
Full documentation, including all API details and examples, is available at:

https://NicolauUP.github.io/MultifractalTools.jl/stable/

📜 License
This package is provided under the MIT License. 


## 📄 Citing

If you use `MultifractalTools.jl` in your research, please cite the following paper:

> Ricardo Oliveira, Nicolau Sobrosa, Pedro Ribeiro, Bruno Amorim, and Eduardo V. Castro (2025). *Local Density of States as a Probe of Multifractality in Quasiperiodic Moiré Materials*. arXiv preprint arXiv:2510.20575.
https://arxiv.org/abs/2510.20575
