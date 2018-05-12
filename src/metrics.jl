SchwarzschildMetric(M) = begin
    scg00(point) = -(1.0 - (2.0M / point[2]))
    scg11(point) = (1.0 - (2.0M / point[2]))^(-1.0)
    scg22(point) = point[2]^(2.0)
    scg33(point) = point[2]^(2.0) * sin(point[3])

    scg00i(point) = -(1.0 - (2.0M / point[2]))^(-1.0)
    scg11i(point) = (1.0 - (2.0M / point[2]))
    scg22i(point) = point[2]^(-2.0)
    scg33i(point) = point[2]^(-2.0) * csc(point[3])

    return Metric([scg00, scg11, scg22, scg33], [scg00i, scg11i, scg22i, scg33i])
end

KerrNewmanMetric(M, a, Q) = begin
    # Metric mappings
    kng00(x) = begin
        r  = x[2]; θ = x[3]
        Δ  = r*r - 2.0*M*r + a*a + Q*Q
        ρ2 = r*r + a*a*cos(θ)*cos(θ)

        return (a*a*sin(θ)*sin(θ) - Δ) / ρ2 
    end

    kng03(x) = begin
        r  = x[2]; θ = x[3]
        Δ  = r*r - 2.0*M*r + a*a + Q*Q
        ρ2 = r*r + a*a*cos(θ)*cos(θ)

        return (sin(θ)*sin(θ) / ρ2)*(a*Δ - a*r*r - a*a*a)
    end

    kng11(x) = begin
        r  = x[2]; θ = x[3]
        Δ  = r*r - 2.0*M*r + a*a + Q*Q
        ρ2 = r*r + a*a*cos(θ)*cos(θ)

        return ρ2 / Δ
    end

    kng22(x) = begin
        r  = x[2]; θ = x[3]
        Δ  = r*r - 2.0*M*r + a*a + Q*Q
        ρ2 = r*r + a*a*cos(θ)*cos(θ)

        return ρ2
    end

    kng33(x) = begin
        r  = x[2]; θ = x[3]
        Δ  = r*r - 2.0*M*r + a*a + Q*Q
        ρ2 = r*r + a*a*cos(θ)*cos(θ)

        return (sin(θ)*sin(θ) / ρ2) * ((a^4 + 2*a*a*r*r + r^4) - a*a*Δ*sin(θ)*sin(θ))
    end

    kngxx(x) = 0.0

    mat = Matrix{Function}((4, 4))
    mat[1, 1] = kng00; mat[1, 2] = kngxx; mat[1, 3] = kngxx; mat[1, 4] = kng03
    mat[2, 1] = kngxx; mat[2, 2] = kng11; mat[2, 3] = kngxx; mat[2, 4] = kngxx
    mat[3, 1] = kngxx; mat[3, 2] = kngxx; mat[3, 3] = kng22; mat[3, 4] = kngxx
    mat[4, 1] = kng03; mat[4, 2] = kngxx; mat[4, 3] = kngxx; mat[4, 4] = kng33

    # Inverse metric mappings
    ikng00(x) = begin
        r   = x[2]; θ = x[3]
        Δ   = r*r - 2.0*M*r + a*a + Q*Q
        ρ2  = r*r + a*a*cos(θ)*cos(θ)
        dem = Δ * (a*a + 2*r*r + a*a*cos(2*θ))^2

        return (-4*ρ2 / dem)*((a*a + r*r)^2 - a*a*Δ*sin(θ)*sin(θ))
    end

    ikng03(x) = begin
        r   = x[2]; θ = x[3]
        Δ   = r*r - 2.0*M*r + a*a + Q*Q
        ρ2  = r*r + a*a*cos(θ)*cos(θ)
        dem = Δ * (a*a + 2*r*r + a*a*cos(2*θ))^2

        return (-4*ρ2 / dem)*(a * (a*a + r*r - Δ))
    end

    ikng11(x) = begin
        r   = x[2]; θ = x[3]
        Δ   = r*r - 2.0*M*r + a*a + Q*Q
        ρ2  = r*r + a*a*cos(θ)*cos(θ)

        return Δ / ρ2
    end

    ikng22(x) = begin
        r   = x[2]; θ = x[3]
        Δ   = r*r - 2.0*M*r + a*a + Q*Q
        ρ2  = r*r + a*a*cos(θ)*cos(θ)

        return 1 / ρ2
    end

    ikng33(x) = begin
        r   = x[2]; θ = x[3]
        Δ   = r*r - 2.0*M*r + a*a + Q*Q
        ρ2  = r*r + a*a*cos(θ)*cos(θ)
        dem = Δ * (a*a - (a*a + r*r)*csc(θ)*csc(θ))^2

        return (ρ2 / dem)*(Δ*csc(θ)*csc(θ) - a*a)*csc(θ)^4
    end

    ikngxx(x) = 0.0
    
    imat = Matrix{Function}((4, 4))
    imat[1, 1] = ikng00; imat[1, 2] = ikngxx; imat[1, 3] = ikngxx; imat[1, 4] = ikng03
    imat[2, 1] = ikngxx; imat[2, 2] = ikng11; imat[2, 3] = ikngxx; imat[2, 4] = ikngxx
    imat[3, 1] = ikngxx; imat[3, 2] = ikngxx; imat[3, 3] = ikng22; imat[3, 4] = ikngxx
    imat[4, 1] = ikng03; imat[4, 2] = ikngxx; imat[4, 3] = ikngxx; imat[4, 4] = ikng33

    return Metric(mat, imat)
end

killing_e(metric::Metric, init_x, init_u) = begin
    gμν, _ = evaluate(metric, init_x)

    return -gμν[1, :]'init_u
end

killing_l(metric::Metric, init_x, init_uϕ) = begin
    gμν, _ = evaluate(metric, init_x)
    K = [0, 0, 0, 1]

    return gμν[4, :]'init_uϕ
end

export killing_e, killing_l
