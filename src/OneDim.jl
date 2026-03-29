"""
    add_exchange1D!(spins, p)

Accumulate nearest-neighbor exchange contributions to `p.fields.Beff` for a 1D
open chain.
"""
function add_exchange1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    @inbounds for i = 1:p.Nx
        if i > 1 # Has neighbor to the left
            for c = 1:3
                p.fields.Beff[c,i] += -p.J * spins[c, i-1]
            end
        end
        if i < p.Nx # Has neighbor to the right
            for c = 1:3
                p.fields.Beff[c,i] += -p.J * spins[c, i+1]
            end
        end
    end
end

"""
    add_bulk_dmi1D!(spins, p)

Accumulate translationally invariant nearest-neighbor bulk DMI contributions
for the discrete bond energy `sum_i D · (S_i × S_{i+1})` on a 1D open chain.
This contributes `D × (S_{i+1} - S_{i-1})` to `p.fields.Beff`.
"""
function add_bulk_dmi1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    Dx, Dy, Dz = p.D

    @inbounds for i = 1:p.Nx
        dSx = (i < p.Nx ? spins[1, i + 1] : 0.0) - (i > 1 ? spins[1, i - 1] : 0.0)
        dSy = (i < p.Nx ? spins[2, i + 1] : 0.0) - (i > 1 ? spins[2, i - 1] : 0.0)
        dSz = (i < p.Nx ? spins[3, i + 1] : 0.0) - (i > 1 ? spins[3, i - 1] : 0.0)

        p.fields.Beff[1, i] += Dy * dSz - Dz * dSy
        p.fields.Beff[2, i] += Dz * dSx - Dx * dSz
        p.fields.Beff[3, i] += Dx * dSy - Dy * dSx
    end
end

"""
    add_anisotropy1D!(spins, p)

Accumulate on-site anisotropy contributions `-K .* S` into `p.fields.Beff`.
"""
function add_anisotropy1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    @inbounds for i = 1:p.Nx
        for c = 1:3
            p.fields.Beff[c,i] += -p.K[c] * spins[c, i]
        end
    end
end

"""
    add_Bext1D!(spins, p)

Accumulate uniform external-field contributions into `p.fields.Beff`.
"""
function add_Bext1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    @inbounds for i = 1:p.Nx
        for c = 1:3
            p.fields.Beff[c,i] += -p.B[c]
        end
    end
end

"""
    add_B_stag1D!(spins, p)

Accumulate a staggered z-field term with alternating sign along the chain.
"""
function add_B_stag1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    @inbounds for i = 1:p.Nx
        sign = isodd(i) ? -1.0 : 1.0
        p.fields.Beff[3,i] += -p.B_stag * sign
    end
end

"""
    add_stt1D!(spins, p)

Accumulate the 1D Zhang-Li spin-transfer torque directly into `p.fields.dS_1`. 
Unlike other functions which accumulate into `p.fields.Beff`, stt goes
directly into the time derivative because the adiabatic term is not of the 
form s x b.

The discrete spatial derivative uses centered differences in the bulk and
one-sided differences at the open boundaries. The added torque is
`-u_stt * ∂x s + beta_stt * u_stt * (s × ∂x s)`.
The implementation assumes lattice-spacing units with `a = 1`.

# References
- Shufeng Zhang and Zhidong Li, "Roles of Nonequilibrium Conduction Electrons
  on the Magnetization Dynamics of Ferromagnets," Phys. Rev. Lett. 93, 127204
  (2004), https://doi.org/10.1103/PhysRevLett.93.127204.
"""
function add_stt1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    if p.Nx == 1 || p.u_stt == 0.0
        return nothing
    end

    prefactor = p.beta_stt * p.u_stt

    @inbounds for i = 1:p.Nx
        if i == 1
            dSx = spins[1, 2] - spins[1, 1]
            dSy = spins[2, 2] - spins[2, 1]
            dSz = spins[3, 2] - spins[3, 1]
        elseif i == p.Nx
            dSx = spins[1, p.Nx] - spins[1, p.Nx - 1]
            dSy = spins[2, p.Nx] - spins[2, p.Nx - 1]
            dSz = spins[3, p.Nx] - spins[3, p.Nx - 1]
        else
            dSx = 0.5 * (spins[1, i + 1] - spins[1, i - 1])
            dSy = 0.5 * (spins[2, i + 1] - spins[2, i - 1])
            dSz = 0.5 * (spins[3, i + 1] - spins[3, i - 1])
        end

        tx = spins[2, i] * dSz - spins[3, i] * dSy
        ty = spins[3, i] * dSx - spins[1, i] * dSz
        tz = spins[1, i] * dSy - spins[2, i] * dSx

        p.fields.dS_1[1, i] += -p.u_stt * dSx + prefactor * tx
        p.fields.dS_1[2, i] += -p.u_stt * dSy + prefactor * ty
        p.fields.dS_1[3, i] += -p.u_stt * dSz + prefactor * tz
    end

    return nothing
end

"""
    add_gilbert1D!(spins, p)

Accumulates local Gilbert damping contributions based on the current `dS_2`.
"""
function add_gilbert1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    @inbounds for i = 1:p.Nx
        for c = 1:3
            p.fields.Beff[c,i] += p.αG * p.fields.dS_2[c,i]
        end
    end
end

"""
    add_nloc_damping1D!(spins, p)

Accumulate translationally invariant nonlocal damping contributions using
`p.ker_dx` and `p.Λtens`.
"""
function add_nloc_damping1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    @inbounds for i = 1:p.Nx
        for (k,dx) in enumerate(p.ker_dx)
            j = i + dx
            if 1 <= j <= p.Nx
                for a=1:3
                    for b=1:3
                        p.fields.Beff[a,i] += p.Λtens[a,b,k] * p.fields.dS_2[b,j]
                    end
                end
            end
        end
    end
#     println("Added nonlocal damping")
end     

"""
    add_nloc_damping_stag1D!(spins, p)

Accumulates staggered nonlocal damping terms using sublattice-dependent tensors.
"""
function add_nloc_damping_stag1D!(spins::Array{Float64, 2}, p::LLGParams1D)
    @inbounds for i = 1:p.Nx
        for (k,dx) in enumerate(p.ker_dx_stag)
            j = i + dx
            if 1 <= j <= p.Nx
                sublat = isodd(i) ? 1 : 2
                for a=1:3
                    for b=1:3
                        p.fields.Beff[a,i] += p.Λtens_stag[sublat,a,b,k] * p.fields.dS_2[b,j]
                    end
                end
            end
        end
    end
end     


"""
    normalize_spins1D!(u, p, t; verbose=false)

Normalize each spin vector in-place to unit length. Intended for use in the
DiscreteCallback during time integration.
"""
function normalize_spins1D!(u, p, t; verbose=false)
    spins = reshape(u, 3, p.Nx)
    normalize_spins!(spins)
    if verbose println("time is $t"); flush(stdout) end
end

"""
    rhs1D!(spins, p, t)

Evaluate the 1D LLG right-hand side for the current spin configuration.

Returns `p.fields.dS_2`, which stores the latest time derivative estimate.
"""
function rhs1D!(spins::Array{Float64, 2}, p::LLGParams1D, t::Float64)
    fill!(p.fields.Beff, 0.0)

    add_exchange1D!(spins, p)
    add_bulk_dmi1D!(spins, p)
    add_anisotropy1D!(spins, p)
    add_Bext1D!(spins, p)
    if p.stag add_B_stag1D!(spins, p) end

    @inbounds for i = 1:p.Nx
        @views cross_inplace!(p.fields.dS_1[:,i], spins[:,i], p.fields.Beff[:,i])
    end
    if p.stt_active add_stt1D!(spins, p) end

    p.fields.dS_2 .= p.fields.dS_1

    for iter in 1:3
        fill!(p.fields.Beff, 0.0) 
        
        add_gilbert1D!(spins, p)
        add_nloc_damping1D!(spins, p)
        if p.stag add_nloc_damping_stag1D!(spins, p) end

        @inbounds for i = 1:p.Nx
            @views cross_inplace!(p.fields.dS_2[:,i], spins[:,i], p.fields.Beff[:,i])
        end

        p.fields.dS_2 .+= p.fields.dS_1
    end

    return p.fields.dS_2
end 

"""
    rhs1D_DE!(du, u, p, t)

DifferentialEquations-compatible RHS wrapper for 1D LLG dynamics.
"""
function rhs1D_DE!(du, u, p, t)
    spins = reshape(u, 3, p.Nx)

    rhs1D!(spins, p, t)

    du .= reshape(p.fields.dS_2, 3*p.Nx)

    return nothing
end
