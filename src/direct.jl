struct DirectGeodesicProblem
    metric::Metric
    Γ::Array{<: Number, 3}
    ∂gμν::Array{<: Number, 3}
    gμν::Matrix{<: Number}
    gμνinv::Matrix{<: Number}
    diff::Vector{<: Number}
    
    function DirectGeodesicProblem(met::Metric, 
                                   Γc::Array{<: Number, 3}, 
                                   ∂gμνc::Array{<: Number, 3}, 
                                   gμνc::Matrix{<: Number}, 
                                   gμνinvc::Matrix{<: Number}, 
                                   diffc::Vector{<: Number})
        new(met, Γc, ∂gμνc, gμνc, gμνinvc, diffc)
    end
end

function DirectGeodesicProblem(met::Metric)
    Γ = zeros(dim(met), dim(met), dim(met))
    ∂gμν = zeros(dim(met), dim(met), dim(met))
    gμν = zeros(size(met))
    gμνinv = zeros(size(met))
    diff = zeros(dim(met))

    return DirectGeodesicProblem(met, Γ, ∂gμν, gμν, gμνinv, diff)
end

function geodesicstep(dX, X, p::DirectGeodesicProblem, τ)
    # xdot^μ = u^μ
    # udot^μ = -Γ[μ, α, β] u^α u^β

    x = @view X[1:4]
    u = @view X[5:8]
    
    dX[1] = u[1]
    dX[2] = u[2]
    dX[3] = u[3]
    dX[4] = u[4]

    # Update the metrics
    evaluate!(p.gμν, p.gμνinv, p.metric, x)

    # Update the connection coefficients at the current point
    christoffel!(p.Γ, p.∂gμν, p.gμν, p.gμνinv, p.diff, p.metric, x)

    dX[5] = -u'p.Γ[:, :, 1]*u
    dX[6] = -u'p.Γ[:, :, 2]*u
    dX[7] = -u'p.Γ[:, :, 3]*u
    dX[8] = -u'p.Γ[:, :, 4]*u
end

export geodesicstep
