struct DirectGeodesicProblem
    metric::Metric
end

function geodesic!(dX, X, p, tau)
    # xdot^μ = u^μ
    # udot^μ = -ch[α, β, μ] u^α u^β

    x = X[1:4]
    u = X[5:8]
    
    dX[1] = u[1]
    dX[2] = u[2]
    dX[3] = u[3]
    dX[4] = u[4]

    # TODO: Make this more efficient with more stuff in the DirectGeodesicProblem later
    ch = christoffel(p.metric, x)

    dX[5] = -u'ch[:, :, 1]*u
    dX[6] = -u'ch[:, :, 2]*u
    dX[7] = -u'ch[:, :, 3]*u
    dX[8] = -u'ch[:, :, 4]*u
end
