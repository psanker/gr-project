__precompile__(true)

module GeodesicIntegration

include("util.jl")
include("metric.jl")
include("adm.jl")

# General functions
export grad, ∂

# Metric related exports
export Metric, 
       dim, 
       evaluate, 
       christoffel

# ADM-specific behaviour
export threemetric, 
       lapse,
       shift, 
       hamiltonian

end
