# Examples

## Notebook files in this repository

- `examples/LambdaLLG_1D_quickstart.ipynb`
- `examples/LambdaLLG_2D_quickstart.ipynb`
- `LambdaLLG_examples.ipynb` (index notebook linking both quickstarts)

## Typical 1D workflow

1. Build `LLGParams1D` with exchange/anisotropy/field/damping settings.
2. Prepare `s0` as a `3 x Nx` matrix of unit spin vectors.
3. Run `evolve1D(s0, tspan, p)`.
4. Convert output with `format_results(sol)` or post-process `sol.u` directly.

## Typical 2D workflow

1. Build `LLGParams2D(Nx, Ny, ...)`.
2. Prepare `s0` as a `3 x Nx x Ny` tensor.
3. Run `evolve2D(s0, tspan, p)`.
4. Analyze observables such as average magnetization versus time.

## Enabling nonlocal damping

Populate kernel offsets and tensors before calling `evolve1D`/`evolve2D`:

```julia
p.ker_dx = collect(-2:2)
p.Λtens = zeros(3, 3, length(p.ker_dx))
```

For 2D, set both `ker_dx`, `ker_dy`, and `Λtens[a,b,kx,ky]`.
