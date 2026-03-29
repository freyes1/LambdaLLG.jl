# LambdaLLG.jl

`LambdaLLG.jl` is a lightweight Julia package for integrating Landau-Lifshitz-Gilbert (LLG)
dynamics with local and nonlocal damping in 1D and 2D spin lattices.

## Features

- 1D and 2D LLG solvers with open boundary conditions.
- Exchange, anisotropy, external magnetic field, 1D bulk DMI, and directional
  2D nearest-neighbor DMI terms.
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

# Optional 1D bulk DMI with a constant bond vector:
# p = LLGParams1D(Nx, 1.0, (0.0, 0.0, -0.02), (0.0, 0.0, -0.1), 0.05; D=(1.0, 0.0, 0.0))

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
- `examples/LambdaLLG_1D_Bloch_domain_wall.ipynb`
- `examples/LambdaLLG_1D_quickstart.ipynb`
- `examples/LambdaLLG_2D_DMI_skyrmion_test.ipynb`
- `examples/LambdaLLG_2D_quickstart.ipynb`
- `examples/LambdaLLG_seed_relaxation_examples.ipynb`

## Uniform seeds

You can build uniformly magnetized starting configurations directly in the
package:

```julia
using LambdaLLG

s1 = uniform_seed1D(101; direction=(0.0, 0.0, 1.0))
s2 = uniform_seed2D(101, 101; direction=(0.0, 0.0, 1.0))

paint_domain_wall1D!(
    s1;
    center=51.0,
    width=4.0,
    domain_direction=(0.0, 0.0, 1.0),
    wall_normal=(0.0, 1.0, 0.0),
)
normalize_spins!(s1)

paint_domain_wall2D!(
    s2;
    point=(51.0, 51.0),
    slope=0.5,
    width=4.0,
    domain_direction=(0.0, 0.0, 1.0),
    wall_normal=(1.0, 0.0, 0.0),
)
normalize_spins!(s2)

paint_skyrmion2D!(
    s2;
    center=(51.0, 51.0),
    radius=12.0,
    width=3.0,
    center_direction=-1,
    helicity=pi / 2,
    vorticity=1.0,
)
normalize_spins!(s2)
```

Use `normalize_spins!` if you modify a seed in-place before passing it to the
solver. The 1D domain-wall painter is additive, so you can paint multiple walls
onto the same uniform background before normalizing. Reverse
`domain_direction` to flip which domain points "up" or "down".
In 2D, `wall_normal = (1, 0, 0)` means local front-normal (Neel-like) and
`wall_normal = (0, 1, 0)` means local along-wall (Bloch-like).
For skyrmions, `helicity = 0` gives a Neel-type texture and `helicity = pi / 2`
gives a Bloch-type texture when `center_direction = -1`, so the
background points along `+z`.

## Notes

- Spins are renormalized at every callback step during integration.
- `format_results(sol)` returns rows as time snapshots.
- For 1D bulk DMI, pass a constant bond vector with the `D=(Dx, Dy, Dz)` keyword.
- For 2D DMI, pass directional bond vectors with
  `D=((Dx_x, Dx_y, Dx_z), (Dy_x, Dy_y, Dy_z))`.
- On a square lattice, common choices are bulk DMI
  `D=((D, 0, 0), (0, D, 0))` and interface DMI
  `D=((0, D, 0), (-D, 0, 0))`.
- For nonlocal damping, populate `ker_dx`, `ker_dy`, and `Λtens` before evolution.
