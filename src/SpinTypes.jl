"""
    LLGFields1D

Workspace buffers for the 1D Landau–Lifshitz–Gilbert solver.

# Fields
- `Beff::Array{Float64, 2}`: Effective magnetic field array of size (3, Nx).
- `dS::Array{Float64, 2}`: Temporary storage for spin time derivatives of size (3, Nx).

These arrays are preallocated and reused during RHS evaluations to avoid allocations.
"""
struct LLGFields1D
    Beff::Array{Float64, 2}
    dS_1::Array{Float64, 2}
    dS_2::Array{Float64, 2}
end

struct LLGFields2D
    Beff::Array{Float64, 3}
    dS_1::Array{Float64, 3}
    dS_2::Array{Float64, 3}
end

"""
    LLGParams1D

Parameters and precomputed data for a 1D LLG simulation with open boundary conditions.

# Fields
- `Nx::Int`: Number of lattice sites.
- `J::Float64`: Nearest-neighbor exchange constant.
- `K::NTuple{3, Float64}`: On-site anisotropy vector.
- `B::NTuple{3, Float64}`: External magnetic field.
- `αG::Float64`: Local Gilbert damping constant.
- `ker_dx::Vector{Int}`: Relative site offsets for nonlocal damping kernel.
- `Λtens::Array{Float64,3}`: Nonlocal damping tensor of size (3,3,n_offsets).
- `fields::LLGFields1D`: Preallocated workspace buffers.

# Notes
- Open boundary conditions are assumed.
- Translation invariance is assumed for the nonlocal kernel.
"""
mutable struct LLGParams1D
    Nx::Int
    J::Float64
    K::NTuple{3, Float64}
    B::NTuple{3, Float64}
    αG::Float64
    ker_dx::Vector{Int}
    Λtens::Array{Float64,3}
    fields::LLGFields1D
    stag::Bool
    B_stag::Float64
    ker_dx_stag::Vector{Int}
    Λtens_stag::Array{Float64, 4}

    function LLGParams1D(
        Nx::Int,
        J::Float64,
        K::NTuple{3, Float64},
        B::NTuple{3, Float64},
        αG::Float64;
        ker_dx::Vector{Int} = Int[],
        Λtens::Array{Float64,3} = zeros(3,3,0)
    )

        Beff = zeros(3, Nx)
        dS_1 = zeros(3, Nx)
        dS_2 = zeros(3, Nx)

        fields = LLGFields1D(Beff, dS_1, dS_2)

        return new(Nx, J, K, B, αG, ker_dx, Λtens, fields, 
                   false, 0.0, Int[], zeros(2,3,3,0))
    end
end

mutable struct LLGParams2D
    Nx::Int
    Ny::Int
    J::Float64
    K::NTuple{3, Float64}
    B::NTuple{3, Float64}
    αG::Float64
    ker_dx::Vector{Int}
    ker_dy::Vector{Int}
    Λtens::Array{Float64,4}
    fields::LLGFields2D
    stag::Bool
    B_stag::Float64
    ker_dx_stag::Vector{Int}
    ker_dy_stag::Vector{Int}
    Λtens_stag::Array{Float64,5}

    function LLGParams2D(
        Nx::Int,
        Ny::Int,
        J::Float64,
        K::NTuple{3, Float64},
        B::NTuple{3, Float64},
        αG::Float64;
        ker_dx::Vector{Int} = Int[],
        ker_dy::Vector{Int} = Int[],
        Λtens::Array{Float64,4} = zeros(3,3,0,0)
    )

        Beff = zeros(3, Nx, Ny)
        dS_1 = zeros(3, Nx, Ny)
        dS_2 = zeros(3, Nx, Ny)

        fields = LLGFields2D(Beff, dS_1, dS_2)

        return new(Nx, Ny, J, K, B, αG, ker_dx, ker_dy, Λtens, fields, 
                   false, 0.0, Int[], Int[], zeros(2,3,3,0,0))
    end
end
