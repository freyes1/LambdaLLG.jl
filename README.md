# LambdaLLG.jl

`LambdaLLG.jl` is a lightweight Julia package for integrating Landau-Lifshitz-Gilbert (LLG)
dynamics with local and nonlocal damping in 1D and 2D spin lattices.

## Features

- 1D and 2D LLG solvers with open boundary conditions.
- Exchange, anisotropy, and external magnetic field terms.
- Local Gilbert damping and optional nonlocal damping kernels.
- Optional staggered (two-sublattice) damping/field terms.
- DifferentialEquations.jl backend with per-step spin normalization callback.

## Installation (local development)

```julia
] activate .
] develop /path/to/LambdaLLG.jl
```

Then in code:

```julia
using LambdaLLG
```

## Quick start (1D)

```julia
using LambdaLLG
using Random
using Statistics

Random.seed!(1)
Nx = 32

p = LLGParams1D(
    Nx,
    1.0,
    (0.0, 0.0, -0.02),
    (0.0, 0.0, -0.1),
    0.05,
)

s0 = zeros(3, Nx)
for i in 1:Nx
    theta = 0.25
    phi = 2pi * rand()
    s0[:, i] .= (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
end

sol = evolve1D(s0, (0.0, 40.0), p)
mz = [mean(reshape(u, 3, Nx)[3, :]) for u in sol.u]

# Convenience output matrix: first column is time
result = format_results(sol)
```

## Documentation

- User docs source: `docs/src/`
- Build docs locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
```

## Example notebooks

- `LambdaLLG_examples.ipynb` (index notebook)
- `examples/LambdaLLG_1D_quickstart.ipynb`
- `examples/LambdaLLG_2D_quickstart.ipynb`

## Notes

- Spins are renormalized at every callback step during integration.
- `format_results(sol)` returns rows as time snapshots.
- For nonlocal damping, populate `ker_dx`, `ker_dy`, and `Λtens` before evolution.
