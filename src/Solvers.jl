"""
    evolve1D(s0, tspan, p::LLGParams1D; reltol=1e-6, abstol=1e-6, verbose=false, stt=false)

Integrates 1D LLG dynamics from initial state `s0` over `tspan`.

Spin normalization is enforced at each accepted step via a
`DiscreteCallback`. Set `stt=true` to include the optional 1D Zhang-Li
spin-transfer torque using the `u_stt` and `beta_stt` fields stored in `p`.
The STT discretization assumes lattice-spacing units with `a = 1`. Returns the
`ODESolution`.
"""
function evolve1D(s0, tspan, p::LLGParams1D; reltol = 1e-6, abstol = 1e-6, verbose=false, stt=false)
    # A callback function could be used to include time dependent parameters
    cb = DE.DiscreteCallback((u,t,integrator)->true,
                      integrator->normalize_spins1D!(integrator.u, integrator.p, integrator.t, verbose=verbose))

    previous_stt_active = p.stt_active
    p.stt_active = stt

    try
        prob = DE.ODEProblem(rhs1D_DE!, reshape(s0, 3*p.Nx), tspan, p)
        return DE.solve(prob, DE.Tsit5(), reltol = reltol, abstol = abstol, callback=cb)
    finally
        p.stt_active = previous_stt_active
    end
end

"""
    evolve2D(s0, tspan, p::LLGParams2D; reltol=1e-6, abstol=1e-6, verbose=false)

Integrate 2D LLG dynamics from initial state `s0` over `tspan`.

Spin normalization is enforced at each accepted step via a
`DiscreteCallback`. Returns the `ODESolution`.
"""
function evolve2D(s0, tspan, p::LLGParams2D; reltol = 1e-6, abstol = 1e-6, verbose=false)
    cb = DE.DiscreteCallback((u,t,integrator)->true,
                      integrator->normalize_spins2D!(integrator.u, integrator.p, integrator.t, verbose=verbose))

    prob = DE.ODEProblem(rhs2D_DE!, reshape(s0, 3*p.Nx*p.Ny), tspan, p)
    sol = DE.solve(prob, DE.Tsit5(), reltol = reltol, abstol = abstol, callback=cb)
    return sol
end
