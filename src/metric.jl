#=
Metric definitions
=#

using TensorOperations
using ForwardDiff

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

function diagmetric!(m1::Matrix{Function}, m2::Matrix{Function}, vec::Vector{Function}, vec2::Vector{Function})
    @assert length(vec) == length(vec2) "The mappings and their inverses must have the same dimensions"
    @assert size(m1) == size(m2) "The allocated matrices must have the same dimensions"
    @assert size(m1)[1] == size(m1)[2] "The allocated matrices must be square"

    for j in 1:length(vec)
        for i in 1:j
            if j == i
                m1[j, i] = vec[i]
            else
                m1[j, i] = _zerof
                m1[i, j] = _zerof
            end
        end
    end

    for j in 1:length(vec)
        for i in 1:j
            if j == i
                m2[j, i] = vec2[i]
            else
                m2[j, i] = _zerof
                m2[i, j] = _zerof
            end
        end
    end

    return Metric(m1, m2)
end

function diagmetric(vec::Vector{Function}, vec2::Vector{Function})
    met    = Matrix{Function}((length(vec), length(vec)))
    invmet = Matrix{Function}((length(vec), length(vec)))

    return diagmetric!(met, invmet, vec, vec2)
end

Metric(vec::Vector{Function}, vec2::Vector{Function}) = diagmetric(vec, vec2)

getindex(metric::Metric, i::I, j::I) where I <: Integer = metric.mappings[i, j]
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

function evaluate(metric::Metric, x::Vector{<: Number})
    evaluated    = zeros(size(metric))
    invevaluated = zeros(size(metric))

    evaluate!(evaluated, invevaluated, metric, x)
    
    return (evaluated, invevaluated)
end

function evaluate!(evaluated::Matrix{<: Number}, invevaluated::Matrix{<: Number}, metric::Metric, x::Vector{<: Number})
    for j in 1:dim(metric)
        for i in 1:j
            if i ≠ j
                evaluated[j, i]    = metric[j, i](x)
                invevaluated[j, i] = inv(metric)[j, i](x)
                
                evaluated[i, j]    = metric[j, i](x)
                invevaluated[i, j] = inv(metric)[j, i](x)
            else
                evaluated[j, i]    = metric[j, i](x)
                invevaluated[j, i] = inv(metric)[j, i](x)
            end
        end
    end
end

function christoffel(metric::Metric, x::Vector{<: Number}) 
    # Evaluate the metric at the point
    gμν, gμνinv = evaluate(metric, x)
    
    # Allocate space for partials
    ∂gμν = zeros(dim(metric), dim(metric), dim(metric))
    Γ    = zeros(dim(metric), dim(metric), dim(metric))
    diff = zeros(length(x))
    
    christoffel!(Γ, ∂gμν, gμν, gμνinv, diff, metric, x)
    return Γ
end

function christoffel!(Γ::Array{<: Number, 3}, 
                      ∂gμν::Array{<: Number, 3},
                      gμν::Matrix{<: Number},
                      gμνinv::Matrix{<: Number},
                      diff::Vector{<: Number},
                      metric::Metric,
                      x::Vector{<: Number})
    
    # Take advantage of symmetries
    for j in 1:dim(metric)
        for i in 1:j
            if i ≠ j
                ForwardDiff.gradient!(diff, metric[j, i], x)
                ∂gμν[j, i, :] = diff
                ∂gμν[i, j, :] = diff
            else
                ForwardDiff.gradient!(diff, metric[j, i], x)
                ∂gμν[j, i, :] = diff
            end
        end
    end
    
    # Split the sum due to limitations of TensorOperations
    @tensoropt begin
        Γ[μ, ν, σ] = 0.5*gμνinv[σ, λ]*∂gμν[λ, μ, ν]
        Γ[μ, ν, σ] = Γ[μ, ν, σ] + 0.5*gμνinv[σ, λ]*∂gμν[λ, ν, μ]
        Γ[μ, ν, σ] = Γ[μ, ν, σ] - 0.5*gμνinv[σ, λ]*∂gμν[μ, ν, λ]
    end
    
    # Clean up numerical weirdness
    cleanartifacts!(Γ)
end
