module LambdaLLG

using LinearAlgebra
using Base.Threads # Currently only for 2D
import DifferentialEquations as DE

export evolve1D, evolve2D
export LLGParams1D, LLGParams2D
export format_results

include("SpinTypes.jl")
include("OneDim.jl")
include("TwoDim.jl")
include("Solvers.jl")

function format_results(sol)
    data = hcat(sol.u...)
    return [sol.t'; data]'
end

function cross_inplace!(out, a, b)
    # Perform manual calculation to avoid allocation
    out[1] = a[2]*b[3] - a[3]*b[2]
    out[2] = a[3]*b[1] - a[1]*b[3]
    out[3] = a[1]*b[2] - a[2]*b[1]
    return nothing
end

end
