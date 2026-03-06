"""
    LambdaLLG

Lightweight Julia package for integrating Landau-Lifshitz-Gilbert (LLG) dynamics
with local and nonlocal damping terms on 1D and 2D lattices.
"""
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

"""
    format_results(sol)

Convert an `ODESolution` into a dense matrix that is convenient for analysis and
plotting.

The returned array has one row per saved time step. The first column is time and
the remaining columns are flattened spin components from the solver state vector.
"""
function format_results(sol)
    data = hcat(sol.u...)
    return [sol.t'; data]'
end

"""
    cross_inplace!(out, a, b)

Compute the cross product `a × b` and write the result into `out` in place.

All vectors are expected to have length 3. This helper is allocation-free and is
used by the RHS kernels.
"""
function cross_inplace!(out, a, b)
    # Perform manual calculation to avoid allocation
    out[1] = a[2]*b[3] - a[3]*b[2]
    out[2] = a[3]*b[1] - a[1]*b[3]
    out[3] = a[1]*b[2] - a[2]*b[1]
    return nothing
end

end
