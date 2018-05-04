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

function sphericaltocartesian(s::Array{<: Number, 1})
    if length(s) == 8
        x = s[2]*sin(s[3])*cos(s[4])
        y = s[2]*sin(s[3])*sin(s[4])
        z = s[2]*cos(s[3])

        ux = s[6]*sin(s[7])*cos(s[8])
        uy = s[6]*sin(s[7])*sin(s[8])
        uz = s[6]*cos(s[7])
        return [s[1], x, y, z, s[5], ux, uy, uz]
    elseif length(s) == 4
        x = s[2]*sin(s[3])*cos(s[4])
        y = s[2]*sin(s[3])*sin(s[4])
        z = s[2]*cos(s[3])
        return [s[1], x, y, z]
    else
        # Assume 3-vec
        x = s[1]*sin(s[2])*cos(s[3])
        y = s[1]*sin(s[2])*sin(s[3])
        z = s[1]*cos(s[2])
        return [x, y, z]
    end
end

function sphericaltocartesian(s::Array{<: Number, 2})
    ret = zeros(size(s))

    for i in 1:size(ret)[1]
        ret[i, :] = sphericaltocartesian(s[i, :])
    end

    return ret 
end

export cleanartifacts!
