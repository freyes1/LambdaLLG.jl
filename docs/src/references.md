# References

This page lists the literature used for the spin-transfer-torque additions and
briefly states what each reference informed.

## Spin-transfer torque form in 1D and 2D

Shufeng Zhang and Zhidong Li, "Roles of Nonequilibrium Conduction Electrons on
the Magnetization Dynamics of Ferromagnets," *Physical Review Letters* **93**,
127204 (2004).

- DOI: <https://doi.org/10.1103/PhysRevLett.93.127204>
- Used for: the phenomenological adiabatic plus nonadiabatic spin-transfer
  torque implemented in this package. In 1D the code uses
  `-u_stt * ∂x m + beta_stt * u_stt * (m × ∂x m)`. In 2D the same local
  Zhang-Li form is implemented with the directional derivative
  `-(u · ∇) m + beta_stt * (m × ((u · ∇) m))`, where `u = (u_x, u_y)`.
  In both cases the discrete derivatives are evaluated in lattice units with
  spacing `a = 1`.

Andre Thiaville, Yuriy Nakatani, Jurgen Miltat, and Nicolas Vernier,
"Micromagnetic understanding of current-driven domain wall motion in patterned
nanowires," *Europhysics Letters* **69**, 990 (2005).

- DOI: <https://doi.org/10.1209/epl/i2004-10452-6>
- Used for: the micromagnetic `u · ∇` rewriting of the local spin-transfer
  torque and the interpretation of the nonadiabatic `beta_stt` term in
  current-driven domain-wall motion.

## Domain-wall sanity check notebook

Andre Thiaville, Yuriy Nakatani, Jurgen Miltat, and Nicolas Vernier,
"Micromagnetic understanding of current-driven domain wall motion in patterned
nanowires," *Europhysics Letters* **69**, 990 (2005).

- DOI: <https://doi.org/10.1209/epl/i2004-10452-6>
- Used for: choosing current-driven domain-wall propagation at zero applied
  field as the qualitative sanity check for the 1D STT notebook.

## Skyrmion-motion sanity check notebook

Junichi Iwasaki, Masahito Mochizuki, and Naoto Nagaosa, "Universal
current-velocity relation of skyrmion motion in chiral magnets," *Nature
Communications* **4**, 1463 (2013).

- DOI: <https://doi.org/10.1038/ncomms2442>
- Used for: motivating current-driven skyrmion translation as the qualitative
  sanity check for the 2D STT notebook and for the expectation that a driven
  skyrmion can show predominantly longitudinal motion together with a smaller
  transverse deflection.
