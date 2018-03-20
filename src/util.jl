#=
Collection of mathematical utility functions
=#

"""
    grad(f::Function, x::Vector{<: Number})

Evaluates the gradient at a point `x` (if `::Vector`)
"""
const testarray = [-1, -1, -1, -1]

function grad(f::Function, point::Vector{<: Number}; h=1e-4, set=testarray)
    partials = zeros(length(point))

    # Raise type of test point to complex
    testpt = complex(copy(point))

    # Check if gradient is over specified set of variables
    if set == testarray
        for j in 1:length(point)
            # Adjust test point
            testpt[j] += h*im
            
            # Evaluate complex step
            partials[j] = imag(f(testpt)) / h

            # Reset test point
            testpt = real(testpt) + 0.0im
        end
    else
        for k in set
            # Adjust test point
            testpt[k] += h*im
            
            # Evaluate complex step
            partials[k] = imag(f(testpt)) / h

            # Reset test point
            testpt = real(testpt) + 0.0im
        end
    end
    
    return partials
end

"""
    grad(f::Function, x::T) where T <: Real

Evaluates the derivative of a univariate function `f` at a point `x` via the complex-step method
"""
grad(f::Function, x::T; h=1e-4) where T <: Real = imag(f(complex(x) + h*im)) / h

"""
    ∂(f::Function, point::Vector{<: Number}; h=1e-4)

Alias to `grad(f, point; h)`
"""
∂(f::Function, point::Vector{<: Number}; h=1e-4, set=testarray) = grad(f, point; h=h, set=set)

"""
    ∂(f::Function, x::T; h=1e-4) where T <: Real

Alias to `grad(f, x; h)`
"""
∂(f::Function, x::T; h=1e-4) where T <: Real = grad(f, x; h=h)

@inbounds function cleanartifacts!(arr::T) where T <: AbstractArray
    for ind in eachindex(arr)
        if abs(arr[ind]) < eps(abs(arr[ind]))
            arr[ind] = zero(arr[ind])
        end
    end
end
