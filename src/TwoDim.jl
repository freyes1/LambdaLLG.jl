"""
    add_exchange2D!(spins, p)

Accumulate nearest-neighbor exchange contributions to `p.fields.Beff` on a 2D
open lattice.
"""
function add_exchange2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny 
        @inbounds for i = 1:p.Nx 
        if i > 1 # Has neighbor to the left
            for c = 1:3
                p.fields.Beff[c,i,j] += -p.J * spins[c, i-1, j]
            end
        end
        if i < p.Nx # Has neighbor to the right
            for c = 1:3
                p.fields.Beff[c,i,j] += -p.J * spins[c, i+1,j]
            end
        end
        if j > 1 # Has neighbor to the below
            for c = 1:3
                p.fields.Beff[c,i,j] += -p.J * spins[c, i,j-1]
            end
        end
        if j < p.Ny # Has neighbor to the above
            for c = 1:3
                p.fields.Beff[c,i,j] += -p.J * spins[c, i,j+1]
            end
        end
    end
    end
end

"""
    add_anisotropy2D!(spins, p)

Accumulate on-site anisotropy contributions `-K .* S` into `p.fields.Beff`.
"""
function add_anisotropy2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            for c = 1:3
                p.fields.Beff[c,i, j] += -p.K[c] * spins[c, i, j]
            end
        end
    end
end

"""
    add_Bext2D!(spins, p)

Accumulate uniform external-field contributions into `p.fields.Beff`.
"""
function add_Bext2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            for c = 1:3
                p.fields.Beff[c,i,j] += -p.B[c]
            end
        end
    end
end

"""
    add_B_stag2D!(spins, p)

Accumulate a checkerboard staggered z-field term with alternating sign between
sublattices.
"""
function add_B_stag2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            sign = isodd(i + j) ? -1.0 : 1.0
            p.fields.Beff[3,i,j] += -p.B_stag * sign
        end
    end
end

"""
    add_gilbert2D!(spins, p)

Accumulate local Gilbert damping contributions based on the current `dS_2`.
"""
function add_gilbert2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            for c = 1:3
                p.fields.Beff[c,i,j] += p.αG * p.fields.dS_2[c,i,j]
            end
        end
    end
end

# function add_nloc_damping2D!(spins::Array{Float64, 3}, p::LLGParams2D)
#     temp = zeros(Float64, 3,1)
#     @inbounds for j = 1:p.Ny
#         @inbounds for i = 1:p.Nx
#         for (ky,dy) in enumerate(p.ker_dy)
#             shiftj = j + dy
#         for (kx,dx) in enumerate(p.ker_dx)
#             shifti = i + dx
            
#             if (1 <= shifti <= p.Nx ) && (1 <= shiftj <= p.Ny)
#                 @views p.fields.Beff[:,i,j] .+= mul!(temp, p.Λtens[:,:,kx,ky], p.fields.dS_2[:,shifti,shiftj])
#             end
#         end
#         end
#     end
#     end
#     println("Added nonlocal damping")
# end  

"""
    add_nloc_damping2D!(spins, p)

Accumulate translationally invariant nonlocal damping contributions using the
offset supports `ker_dx` and `ker_dy`.
"""
function add_nloc_damping2D!(spins::Array{Float64,3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            for (ky,dy) in enumerate(p.ker_dy)
                shiftj = j + dy
                if !(1 ≤ shiftj ≤ p.Ny); continue; end

                for (kx,dx) in enumerate(p.ker_dx)
                    shifti = i + dx
                    if !(1 ≤ shifti ≤ p.Nx); continue; end
                    
                    for a = 1:3
                        for b = 1:3
                            p.fields.Beff[a,i,j] +=  p.Λtens[a,b,kx,ky]*p.fields.dS_2[b,shifti,shiftj]
                        end
                    end
                end
            end
        end
    end
#     println("Added nonlocal damping")
end

"""
    add_nloc_damping_stag2D!(spins, p)

Accumulate sublattice-dependent nonlocal damping contributions in 2D.
"""
function add_nloc_damping_stag2D!(spins::Array{Float64,3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            for (ky,dy) in enumerate(p.ker_dy_stag)
                shiftj = j + dy
                if !(1 ≤ shiftj ≤ p.Ny); continue; end

                for (kx,dx) in enumerate(p.ker_dx_stag)
                    shifti = i + dx
                    if !(1 ≤ shifti ≤ p.Nx); continue; end
                    
                    sublat = isodd(i+j) ? 1 : 2

                    for a = 1:3
                        for b = 1:3
                            p.fields.Beff[a,i,j] +=  p.Λtens_stag[sublat,a,b,kx,ky]*p.fields.dS_2[b,shifti,shiftj]
                        end
                    end
                end
            end
        end
    end
#     println("Added nonlocal damping")
end



"""
    normalize_spins2D!(u, p, t; verbose=false)

Normalize each 2D lattice spin vector in-place to unit length. Intended for use
in the DiscreteCallback during time integration.
"""
function normalize_spins2D!(u, p, t; verbose=false)
    spins = reshape(u, 3, p.Nx, p.Ny)

    @threads for j in 1:p.Ny
        @inbounds for i in 1:p.Nx
            n = sqrt(spins[1,i,j]^2 + spins[2,i,j]^2 + spins[3,i,j]^2)
            for c in 1:3
                spins[c,i,j] /= n
            end
        end
    end
    
    if verbose println("time is $t"); flush(stdout) end
end
    
"""
    rhs2D!(spins, p, t)

Evaluate the 2D LLG right-hand side for the current spin configuration.

Returns `p.fields.dS_2`, which stores the latest time derivative estimate.
"""
function rhs2D!(spins::Array{Float64, 3}, p::LLGParams2D, t::Float64)
    fill!(p.fields.Beff, 0.0)

    add_exchange2D!(spins, p)
    add_anisotropy2D!(spins, p)
    add_Bext2D!(spins, p)
    if p.stag add_B_stag2D!(spins, p) end
    
    @threads for j = 1:p.Ny    
        @inbounds for i = 1:p.Nx
            @views cross_inplace!(p.fields.dS_1[:,i,j], spins[:,i,j], p.fields.Beff[:,i,j])
        end
    end

    p.fields.dS_2 .= p.fields.dS_1

    for iter in 1:3
        fill!(p.fields.Beff, 0.0) 
        
        add_gilbert2D!(spins, p)
        add_nloc_damping2D!(spins, p)
        if p.stag add_nloc_damping_stag2D!(spins, p) end

        @threads for j = 1:p.Ny
            @inbounds for i = 1:p.Nx
                @views cross_inplace!(p.fields.dS_2[:,i,j], spins[:,i,j], p.fields.Beff[:,i,j])
            end
        end

        p.fields.dS_2 .+= p.fields.dS_1
    end

    return p.fields.dS_2
end 


"""
    rhs2D_DE!(du, u, p, t)

DifferentialEquations-compatible RHS wrapper for 2D LLG dynamics.
"""
function rhs2D_DE!(du, u, p, t)
    spins = reshape(u, 3, p.Nx, p.Ny)

    rhs2D!(spins, p, t)

    du .= reshape(p.fields.dS_2, 3*p.Nx*p.Ny)

    return du
end
