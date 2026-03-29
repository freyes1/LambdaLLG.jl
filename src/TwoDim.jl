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
    add_dmi2D!(spins, p)

Accumulate translationally invariant nearest-neighbor DMI contributions for the
discrete bond energy

`sum_{i,j} D_x · (S[i,j] × S[i+1,j]) + D_y · (S[i,j] × S[i,j+1])`

on a 2D open lattice. The two bond vectors in `p.D` correspond to the x- and
y-directed nearest-neighbor bonds, so this same representation covers both bulk
and interface DMI.
"""
function add_dmi2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    Dx, Dy = p.D
    Dxx, Dxy, Dxz = Dx
    Dyx, Dyy, Dyz = Dy

    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            dSx_x = (i < p.Nx ? spins[1, i + 1, j] : 0.0) - (i > 1 ? spins[1, i - 1, j] : 0.0)
            dSy_x = (i < p.Nx ? spins[2, i + 1, j] : 0.0) - (i > 1 ? spins[2, i - 1, j] : 0.0)
            dSz_x = (i < p.Nx ? spins[3, i + 1, j] : 0.0) - (i > 1 ? spins[3, i - 1, j] : 0.0)

            dSx_y = (j < p.Ny ? spins[1, i, j + 1] : 0.0) - (j > 1 ? spins[1, i, j - 1] : 0.0)
            dSy_y = (j < p.Ny ? spins[2, i, j + 1] : 0.0) - (j > 1 ? spins[2, i, j - 1] : 0.0)
            dSz_y = (j < p.Ny ? spins[3, i, j + 1] : 0.0) - (j > 1 ? spins[3, i, j - 1] : 0.0)

            p.fields.Beff[1, i, j] += Dxy * dSz_x - Dxz * dSy_x + Dyy * dSz_y - Dyz * dSy_y
            p.fields.Beff[2, i, j] += Dxz * dSx_x - Dxx * dSz_x + Dyz * dSx_y - Dyx * dSz_y
            p.fields.Beff[3, i, j] += Dxx * dSy_x - Dxy * dSx_x + Dyx * dSy_y - Dyy * dSx_y
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
    add_stt2D!(spins, p)

Accumulate the optional 2D Zhang-Li spin-transfer torque directly into
`p.fields.dS_1`.

Unlike the field-based terms accumulated into `p.fields.Beff`, the spin-transfer
torque contributes directly to the time derivative. In lattice-spacing units
with `a = 1`, the implemented 2D form is

`-(u · ∇) s + beta_stt * (s × ((u · ∇) s))`

with `u = (u_x, u_y) = p.u_stt`. Centered finite differences are used in the
bulk, and one-sided differences are used at open boundaries.

References:
- Shufeng Zhang and Zhidong Li, Phys. Rev. Lett. 93, 127204 (2004),
  doi:10.1103/PhysRevLett.93.127204
- Andre Thiaville, Yuriy Nakatani, Jurgen Miltat, and Nicolas Vernier,
  Europhys. Lett. 69, 990 (2005), doi:10.1209/epl/i2004-10452-6
"""
function add_stt2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    ux, uy = p.u_stt
    if (ux == 0.0 && uy == 0.0) || (p.Nx == 1 && p.Ny == 1)
        return nothing
    end

    beta = p.beta_stt

    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            dSx_dx = if p.Nx == 1
                0.0
            elseif i == 1
                spins[1, 2, j] - spins[1, 1, j]
            elseif i == p.Nx
                spins[1, p.Nx, j] - spins[1, p.Nx - 1, j]
            else
                0.5 * (spins[1, i + 1, j] - spins[1, i - 1, j])
            end
            dSy_dx = if p.Nx == 1
                0.0
            elseif i == 1
                spins[2, 2, j] - spins[2, 1, j]
            elseif i == p.Nx
                spins[2, p.Nx, j] - spins[2, p.Nx - 1, j]
            else
                0.5 * (spins[2, i + 1, j] - spins[2, i - 1, j])
            end
            dSz_dx = if p.Nx == 1
                0.0
            elseif i == 1
                spins[3, 2, j] - spins[3, 1, j]
            elseif i == p.Nx
                spins[3, p.Nx, j] - spins[3, p.Nx - 1, j]
            else
                0.5 * (spins[3, i + 1, j] - spins[3, i - 1, j])
            end

            dSx_dy = if p.Ny == 1
                0.0
            elseif j == 1
                spins[1, i, 2] - spins[1, i, 1]
            elseif j == p.Ny
                spins[1, i, p.Ny] - spins[1, i, p.Ny - 1]
            else
                0.5 * (spins[1, i, j + 1] - spins[1, i, j - 1])
            end
            dSy_dy = if p.Ny == 1
                0.0
            elseif j == 1
                spins[2, i, 2] - spins[2, i, 1]
            elseif j == p.Ny
                spins[2, i, p.Ny] - spins[2, i, p.Ny - 1]
            else
                0.5 * (spins[2, i, j + 1] - spins[2, i, j - 1])
            end
            dSz_dy = if p.Ny == 1
                0.0
            elseif j == 1
                spins[3, i, 2] - spins[3, i, 1]
            elseif j == p.Ny
                spins[3, i, p.Ny] - spins[3, i, p.Ny - 1]
            else
                0.5 * (spins[3, i, j + 1] - spins[3, i, j - 1])
            end

            gx = ux * dSx_dx + uy * dSx_dy
            gy = ux * dSy_dx + uy * dSy_dy
            gz = ux * dSz_dx + uy * dSz_dy

            sx = spins[1, i, j]
            sy = spins[2, i, j]
            sz = spins[3, i, j]

            tx = sy * gz - sz * gy
            ty = sz * gx - sx * gz
            tz = sx * gy - sy * gx

            p.fields.dS_1[1, i, j] += -gx + beta * tx
            p.fields.dS_1[2, i, j] += -gy + beta * ty
            p.fields.dS_1[3, i, j] += -gz + beta * tz
        end
    end

    return nothing
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
    normalize_spins!(spins)
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
    add_dmi2D!(spins, p)
    add_anisotropy2D!(spins, p)
    add_Bext2D!(spins, p)
    if p.stag add_B_stag2D!(spins, p) end
    
    @threads for j = 1:p.Ny    
        @inbounds for i = 1:p.Nx
            @views cross_inplace!(p.fields.dS_1[:,i,j], spins[:,i,j], p.fields.Beff[:,i,j])
        end
    end
    if p.stt_active add_stt2D!(spins, p) end

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
