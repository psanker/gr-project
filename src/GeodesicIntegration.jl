__precompile__(true)

module GeodesicIntegration

include("util.jl")
include("metric.jl")
include("adm.jl")
include("direct.jl")

# General functions

# Metric related exports
export Metric, 
       dim, 
       evaluate, 
       christoffel,
       christoffel!

# ADM-specific behaviour
export threemetric, 
       lapse,
       shift, 
       hamiltonian

# Direct integration with known analytic metric from 4D manifold
export DirectGeodesicProblem,
       geodesic!

end
