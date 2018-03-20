#=
Functions directly dealing with a general metric
=#

using TensorOperations

import Base: size, getindex, ndims, inv

_zerof(x...) = 0.0

struct Metric
    mappings::Matrix{Function}
    inversemappings::Matrix{Function}

    function Metric(mat::Matrix{Function}, mat2::Matrix{Function})
        @assert size(mat)[1] == size(mat)[2] "Metrics are square"
        @assert size(mat2)[1] == size(mat2)[2] "Metrics are square"
        @assert size(mat) == size(mat2) "The mappings and their inverses must have the same dimensions"
        
        # I'm going to assume that the mappings are symmetric for now I guess
        new(mat, mat2)
    end
end

function diagmetric(vec::Vector{Function}, vec2::Vector{Function})
    @assert length(vec) == length(vec2) "The mappings and their inverses must have the same dimensions"
    
    met    = Matrix{Function}((length(vec), length(vec)))
    invmet = Matrix{Function}((length(vec), length(vec)))
    
    for j in 1:length(vec), i in 1:length(vec)
        if j == i
            met[j, i] = vec[i]
        else
            met[j, i] = _zerof
        end
    end

    for j in 1:length(vec), i in 1:length(vec)
        if j == i
            invmet[j, i] = vec2[i]
        else
            invmet[j, i] = _zerof
        end
    end
    
    return Metric(met, invmet)
end

Metric(vec::Vector{Function}, vec2::Vector{Function}) = diagmetric(vec, vec2)

getindex(metric::Metric, i::Integer, j::Integer) = metric.mappings[i, j]
getindex(metric::Metric, r1::T, r2::T) where T <: AbstractVector = getindex(metric.mappings, r1, r2)

size(metric::Metric) = size(metric.mappings)
size(metric::Metric, i::Integer) = size(metric.mappings, i)
ndims(metric::Metric) = ndims(metric.mappings)

inv(metric::Metric) = metric.inversemappings

dim(metric::Metric) = size(metric)[1]
function dim(mat::Matrix)
    @assert size(mat)[1] == size(mat)[2] "dim(::Matrix) is only valid for square matrices"
    
    return size(mat)[1]
end

function evaluate(metric::Metric, point::Vector{<: Number})
    evaluated    = zeros(size(metric))
    invevaluated = zeros(size(metric))
    
    for j in 1:dim(metric), i in 1:dim(metric)
        evaluated[j, i]    = metric[j, i](point)
        invevaluated[j, i] = inv(metric)[j, i](point)
    end
    
    return (evaluated, invevaluated)
end

function ∂(f::Function, metric::Metric, point::Vector{<: Number})
    # Allocate space for partials
    
end

function christoffel(metric::Metric, point::Vector{<: Number}) 
    # Evaluate the metric at the point
    gdμν, _ = evaluate(metric, point)
    
    # Allocate space for partials
    ∂gμν = zeros(dim(metric), dim(metric), dim(metric))
    
    for j in 1:dim(metric), i in 1:dim(metric)
        ∂gμν[j, i, :] = ∂(metric[j, i], point)
    end
    
    gμνinv = inv(gμν)
    
    # Split the sum due to limitations of TensorOperations
    @tensoropt begin
        Γ[σ, μ, ν] := 0.5*gμνinv[σ, λ]*∂gμν[λ, μ, ν]
        Γ[σ, μ, ν] = Γ[σ, μ, ν] + 0.5*gμνinv[σ, λ]*∂gμν[λ, ν, μ]
        Γ[σ, μ, ν] = Γ[σ, μ, ν] - 0.5*gμνinv[σ, λ]*∂gμν[μ, ν, λ]
    end
    
    # Clean up numerical weirdness
    cleanartifacts!(Γ)
    return Γ
end
