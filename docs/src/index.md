# MultifractalTools.jl

Welcome to the documentation for `MultifractalTools.jl`.

This package provides a set of tools for performing multifractal analysis.

## Installation

```julia
using Pkg
Pkg.add("MultifractalTools")
```

## Example Usage

Here is a quick example of how to use the package. We can also embed plots! Here we show the "blind" workflow where you don't use the interactive window to define the fitting region!

```@example
using MultifractalTools
using GLMakie

# Create some sample data
data = ones(1024, 1024)

# 1. Get q values
qs = obtain_qs(-2, 2, 41)

# 2. Compute scaling quantities
scaling_data = compute_scaling_quantities(data, qs)

# 3. Compute the spectrum (using a subset of scaling)
spec = compute_spectrum(scaling_data, qs, 4, 20)

# 4. Plot the spectrum
# We call the plotting function directly
plot_spectrum(spec, which=:Spectrum)
```
Here we show the second workflow, using the interactive approach!
```@example 
using MultifractalTools
using GLMakie
data = ones(1024,1024)

# 1. Get q values
qs = obtain_qs(-2, 2, 51)

# 2. Compute scaling quantities
scaling_data = compute_scaling_quantities(data, qs)

# 3. PLot the scaling quantities to analyze the fitting interactively!
fig = plot_to_fit(scaling_data, qs)
display(fig) #This should open a new window if one is on a jupyter notebook!

# 4. Compute the spectrum (using a subset of scaling)
spec = compute_spectrum(scaling_data, qs, 4, 20)

# 4. Plot the spectrum
plot_spectrum(spec, which=:Both)

```
## Page Index

```@index
```