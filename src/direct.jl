struct DirectGeodesicProblem
    metric::Metric
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

    # TODO: Make this more efficient with more stuff in the DirectGeodesicProblem later
    Γ = christoffel(p.metric, x)

    dX[5] = -u'Γ[:, :, 1]*u
    dX[6] = -u'Γ[:, :, 2]*u
    dX[7] = -u'Γ[:, :, 3]*u
    dX[8] = -u'Γ[:, :, 4]*u
end

export geodesicstep
