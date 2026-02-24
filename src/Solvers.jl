function evolve1D(s0, tspan, p::LLGParams1D; reltol = 1e-6, abstol = 1e-6, verbose=false)
    # A callback function could be used to include time dependent parameters
    cb = DE.DiscreteCallback((u,t,integrator)->true,
                      integrator->normalize_spins1D!(integrator.u, integrator.p, integrator.t, verbose))

    prob = DE.ODEProblem(rhs1D_DE!, reshape(s0, 3*p.Nx), tspan, p)
    sol = DE.solve(prob, DE.Tsit5(), reltol = reltol, abstol = abstol, callback=cb)
    return sol
end

function evolve2D(s0, tspan, p::LLGParams2D; reltol = 1e-6, abstol = 1e-6, verbose=false)
    cb = DE.DiscreteCallback((u,t,integrator)->true,
                      integrator->normalize_spins2D!(integrator.u, integrator.p, integrator.t, verbose))

    prob = DE.ODEProblem(rhs2D_DE!, reshape(s0, 3*p.Nx*p.Ny), tspan, p)
    sol = DE.solve(prob, DE.Tsit5(), reltol = reltol, abstol = abstol, callback=cb)
    return sol
end
