# LambdaLLG.jl

`LambdaLLG.jl` integrates Landau-Lifshitz-Gilbert (LLG) spin dynamics in 1D and 2D,
including local Gilbert damping and optional nonlocal damping kernels.

## Features

- 1D and 2D lattice LLG evolution with open boundaries.
- Local and nonlocal damping terms.
- Optional staggered (two-sublattice) fields and damping kernels.
- DifferentialEquations.jl integration with per-step spin normalization.

## Installation

For local development:

```julia
] activate .
] develop /path/to/LambdaLLG.jl
```

Then:

```julia
using LambdaLLG
```

## Quick Start (1D)

```julia
using LambdaLLG
using Random

Random.seed!(4)
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
data = format_results(sol)
```

`data[:, 1]` contains time. The remaining columns store flattened spin components
at each saved time step.

## Output Shape

- 1D state vectors are stored as length `3*Nx`.
- 2D state vectors are stored as length `3*Nx*Ny`.
- `format_results(sol)` returns a matrix with shape `(nt, 1 + nstate)`.

## Next Steps

- Use [`examples.md`](examples.md) for notebook-oriented workflows.
- Use [`api.md`](api.md) for exported API references.
