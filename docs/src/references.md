# References

This page lists the literature used for the 1D spin-transfer-torque addition and
briefly states what each reference informed.

## Spin-transfer torque form

Shufeng Zhang and Zhidong Li, "Roles of Nonequilibrium Conduction Electrons on
the Magnetization Dynamics of Ferromagnets," *Physical Review Letters* **93**,
127204 (2004).

- DOI: <https://doi.org/10.1103/PhysRevLett.93.127204>
- Used for: the phenomenological adiabatic plus nonadiabatic spin-transfer
  torque form implemented in 1D,
  `-u_stt * ∂x m + beta_stt * u_stt * (m × ∂x m)`.
  In this package the discrete derivative is evaluated in lattice units with
  spacing `a = 1`.

## Domain-wall sanity check interpretation

André Thiaville, Yuriy Nakatani, Jürgen Miltat, and Nicolas Vernier,
"Micromagnetic understanding of current-driven domain wall motion in patterned
nanowires," *Europhysics Letters* **69**, 990 (2005).

- DOI: <https://doi.org/10.1209/epl/i2004-10452-6>
- Used for: choosing the current-driven domain-wall setting as the qualitative
  sanity check and for interpreting the role of the nonadiabatic `beta_stt`
  term in steady wall propagation.
