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
    add_gilbert1D!(spins, p)

Accumulate local Gilbert damping contributions based on the current `dS_2`.
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

Accumulate staggered nonlocal damping terms using sublattice-dependent tensors.
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

    @inbounds for i in 1:p.Nx
        n = sqrt(spins[1,i]^2 + spins[2,i]^2 + spins[3,i]^2)
        for c in 1:3
            spins[c,i] /= n
        end
    end
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
    add_anisotropy1D!(spins, p)
    add_Bext1D!(spins, p)
    if p.stag add_B_stag1D!(spins, p) end

    @inbounds for i = 1:p.Nx
        @views cross_inplace!(p.fields.dS_1[:,i], spins[:,i], p.fields.Beff[:,i])
    end

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
