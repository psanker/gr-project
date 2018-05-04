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

end
