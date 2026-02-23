# LambdaLLG

A lightweight Julia package for simulating the Landau–Lifshitz–Gilbert (LLG) equation with exchange, anisotropy, external fields, and nonlocal damping in 1D and 2D spin systems.

---

## Installation (Local Development)

Clone the repository:

```bash
git clone https://github.com/yourusername/LambdaLLG.git
```

To use in a project, run from within folder where you want to run:

```Julia
] activate .
] develop /path/to/LambdaLLG
```

Then from within notebook or code

```Julia
using LambdaLLG
```

## Notes

- Open boundary conditions are currently implemented.

- Spin normalization is enforced via callback function at every time step.

- Designed for grids up to ~100×100.


