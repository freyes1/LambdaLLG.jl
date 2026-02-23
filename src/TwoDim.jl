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

function add_anisotropy2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            for c = 1:3
                p.fields.Beff[c,i, j] += -p.K[c] * spins[c, i, j]
            end
        end
    end
end

function add_Bext2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            for c = 1:3
                p.fields.Beff[c,i,j] += -p.B[c]
            end
        end
    end
end

function add_Bstag2D!(spins::Array{Float64, 3}, p::LLGParams2D)
    @threads for j = 1:p.Ny
        @inbounds for i = 1:p.Nx
            sign = isodd(i + j) ? -1.0 : 1.0
            p.fields.Beff[3,i,j] += -p.Bstag * sign
        end
    end
end

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


function normalize_spins2D!(u, p, t)
    spins = reshape(u, 3, p.Nx, p.Ny)

    @threads for j in 1:p.Ny
        @inbounds for i in 1:p.Nx
            n = sqrt(spins[1,i,j]^2 + spins[2,i,j]^2 + spins[3,i,j]^2)
            for c in 1:3
                spins[c,i,j] /= n
            end
        end
    end
    
    println("time is $t"); flush(stdout)
end
    
function rhs2D!(spins::Array{Float64, 3}, p::LLGParams2D, t::Float64)
    fill!(p.fields.Beff, 0.0)

    add_exchange2D!(spins, p)
    add_anisotropy2D!(spins, p)
    add_Bext2D!(spins, p)
    if p.stag add_Bstag2D!(spins, p) end
    
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

        @threads for j = 1:p.Ny
            @inbounds for i = 1:p.Nx
                @views cross_inplace!(p.fields.dS_2[:,i,j], spins[:,i,j], p.fields.Beff[:,i,j])
            end
        end

        p.fields.dS_2 .+= p.fields.dS_1
    end

    return p.fields.dS_2
end 


function rhs2D_DE!(du, u, p, t)
    spins = reshape(u, 3, p.Nx, p.Ny)

    rhs2D!(spins, p, t)

    du .= reshape(p.fields.dS_2, 3*p.Nx*p.Ny)

    return du
end
