__precompile__(true)

module GeodesicIntegration

    using Reexport
    @reexport using DifferentialEquations # So we don't have to include 'using DifferentialEquations'

    include("util.jl")
    include("metric.jl")
    include("adm.jl")
    include("direct.jl")
    include("metrics.jl")

    # General functions
    export sphericaltocartesian

    # Metric related exports
    export Metric

    # Direct integration with known analytic metric from 4D manifold
    export DirectGeodesicProblem

    function scaffold(problem::DirectGeodesicProblem, initx::Vector{<: Number}, init3u::Vector{<: Number}, λspan; timelike=true)
        @assert length(init3u) == 3 "Give an initial 3-velocity; the time component will be determined by normalization of the 4-vel"

        ϵ      = Int(timelike)                   # If null (photons), will be 0
        gμν, _ = evaluate(problem.metric, initx) # Needed to calculate initial 4-vel
        gij    = @view gμν[2:end, 2:end]

        b  = (gμν[1, 2:end]'init3u)/gμν[1, 1]
        ac = (init3u'gij*init3u + ϵ)/gμν[1, 1]

        u0  = -b + sqrt(b*b - ac)
        X0  = [initx..., u0, init3u...]

        return ODEProblem(geodesicstep, X0, λspan, problem), X0
    end

    # For user codes
    export scaffold

    # Some sample metrics
    export SchwarzschildMetric, KerrNewmanMetric

    function __init__()
        # Help speed up render times by running some user code so __precompile__ can catch it

        M = 3.0
        scmetric = SchwarzschildMetric(M)

        problem = DirectGeodesicProblem(scmetric)
        init_x  = [0.0, 6.0M, π/2, π/4]; # Initial x^μ
        init_u3 = [-1.0, 0.003, 0.08];   # Initial u^i, u^0 will be determined by normalization

        tspan      = (0.0, 5.0)
        diffeq, X0 = scaffold(problem, init_x, init_u3, tspan)

        ε = killing_e(scmetric, init_x, X0[5:8])
        l = killing_l(scmetric, init_x, X0[5:8])

        # Assert energy and angular momentum conservation via manifold projection callback
        # Maintain the normalization of the 4-velocity
        function g(r, u, p, t)
            r[1] = 0; r[2] = 0; r[3] = 0; r[4] = 0;
            r[5] = @. ((1 - (2M / u[2])) * u[5]) - ε
            r[6] = 0; r[7] = 0;
            r[8] = @. u[2]*u[2]*sin(u[3])*sin(u[3])*u[8] - l
        end

        cb  = ManifoldProjection(g) 
        sol = solve(diffeq, reltol=1e-6, callback=cb);
    end
end
