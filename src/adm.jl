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
Shift vector of the 4-metric
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
