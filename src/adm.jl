#=
Functions related to the ADM formalism directly
=#

"""
Extracts the spacelike 3-metric from the spacetime 4-metric
"""
function threemetric(metric::Metric)
    @assert(dim(metric) == 4, "Extracting the 3-metric is only valid for a 4-metric")
    
    Metric(metric[2:end, 2:end], inv(metric)[2:end, 2:end])
end

function threemetric(gμν::Matrix{Float64})
    @assert(dim(gμν) == 4, "Extracting the 3-metric is only valid for a 4-metric")
    
    return gμν[2:end, 2:end]
end

"""
Lapse function of the 4-metric
"""
function lapse(metric::Metric, point::Vector{<: Number})
    # Assert that the dimensionality of the metric is 4
    @assert(dim(metric) == 4, "This method is only defined for the 4-metric")
    
    # Evalute the metric at the point
    gμν, _ = evaluate(metric, point)

    return (-inv(gμν)[1, 1])^(-1//2)
end

function lapse(gμν::Matrix{Float64})
    @assert(dim(gμν) == 4, "This method is only defined for the 4-metric")
    
    return (-inv(gμν)[1, 1])^(-1//2)
end

"""
Shift 3-vector of the 4-metric
"""
function shift(metric::Metric, point::Vector{<: Number})
    @assert(dim(metric) == 4, "This method is only defined for the 4-metric")
    
    # Evalute the metrics at the point
    gμν, _ = evaluate(metric, point)
    gij, _ = evaluate(threemetric(metric), point)
    
    β_j = @view gμν[:, 1][2:end] # No unnecessary copying
    
    return inv(gij) * β_j
end

function shift(gμν::Matrix{Float64})
    @assert(dim(gμν) == 4, "This method is only defined for the 4-metric")
    
    β_j = @view gμν[:, 1][2:end] # No unnecessary copying
    return inv(threemetric(gμν)) * β_j
end

"""
    hamiltonian(metric::Metric, x::Vector{<: Number}, u::Vector{<: Number}; ε=1.0)

Evaluates the Hamiltonian at a specific point `x` (relative to a *coordinate basis*) and a velocity *one-form* `u`

- `metric` : The spacetime 4-metric
- `x` : The position vector
- `u` : The velocity *one-form*
- `ε` : `true` if the test particle is *time-like*; `false` if the test particle is *null-like*
"""
function hamiltonian(metric::Metric, x::Vector{<: Number}, u::Vector{<: Number}; ε=true)
    @assert length(x) == 4 "Position is a 4-vector"
    @assert length(u) == 4 "Velocity is a 4-dual"
    @assert dim(metric) == 4 "This method is only defined for the 4-metric"

    u3    = @view u[2:end] # For summing ease

    γjk, invγjk = evaluate(threemetric(metric), x)
    
    return lapse(metric, x)*(u3'invγjk*u3 + Float64(ε))^(1//2) + shift(metric, x)'u3
end

function hamiltonian(gμν::Matrix{Float64}, x::Vector{<: Number}, u::Vector{<: Number}; ε=true)
    @assert length(x) == 4 "Position is a 4-vector"
    @assert length(u) == 4 "Velocity is a 4-dual"
    @assert dim(gμν) == 4 "This method is only defined for the 4-metric"

    u3 = @view u[2:end] # For convenience

    return lapse(gμν)*(u3'inv(threemetric(gμν))*u3 + Float64(ε))^(1//2) + shift(gμν)'u3
end
