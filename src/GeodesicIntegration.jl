__precompile__(true)

module GeodesicIntegration

    using Reexport
    @reexport using DifferentialEquations # So we don't have to include 'using DifferentialEquations'

    include("util.jl")
    include("metric.jl")
    include("adm.jl")
    include("direct.jl")

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

        u0  = sqrt((-ϵ - init3u'gμν[2:end, 2:end]*init3u) / gμν[1, 1])

        X0 = [initx..., u0, init3u...] # Vector of initial conditions
        return ODEProblem(geodesicstep, X0, λspan, problem)
    end

    # For user codes
    export scaffold


    function __init__()
        # Help speed up render times by running some user code so __precompile__ can catch it

        M = 3.0

        scg00(point) = -(1.0 - (2.0M / point[2]))
        scg11(point) = (1.0 - (2.0M / point[2]))^(-1.0)
        scg22(point) = point[2]^(2.0)
        scg33(point) = point[2]^(2.0) * sin(point[3])

        scg00i(point) = -(1.0 - (2.0M / point[2]))^(-1.0)
        scg11i(point) = (1.0 - (2.0M / point[2]))
        scg22i(point) = point[2]^(-2.0)
        scg33i(point) = point[2]^(-2.0) * csc(point[3])

        scmetric = Metric([scg00, scg11, scg22, scg33], [scg00i, scg11i, scg22i, scg33i]);


        problem = DirectGeodesicProblem(scmetric)
        init_x  = [0.0, 6.0M, π/2, π/4]; # Initial x^μ
        init_u3 = [-1.0, 0.003, 0.08];   # Initial u^i, u^0 will be determined by normalization

        tspan  = (0.0, 5.0)
        diffeq = scaffold(problem, init_x, init_u3, tspan)

        sol = solve(diffeq, reltol=1e-6);
    end

end
