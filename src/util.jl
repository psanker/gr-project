#=
Collection of mathematical utility functions
=#

@inbounds function cleanartifacts!(arr::T) where T <: AbstractArray
    for ind in eachindex(arr)
        if abs(arr[ind]) < eps(abs(arr[ind]))
            arr[ind] = zero(arr[ind])
        end
    end
end
