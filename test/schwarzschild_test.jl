using GeodesicIntegration, Base.Test

# Init
M = 3.0

scg00(point::Vector{<: Number}) = -(1.0 - (2.0M / point[2]))
scg11(point::Vector{<: Number}) = (1.0 - (2.0M / point[2]))^(-1.0)
scg22(point::Vector{<: Number}) = point[2]^(2.0)
scg33(point::Vector{<: Number}) = point[2]^(2.0) * sin(point[3])

scg00i(point::Vector{<: Number}) = -(1.0 - (2.0M / point[2]))^(-1.0)
scg11i(point::Vector{<: Number}) = (1.0 - (2.0M / point[2]))
scg22i(point::Vector{<: Number}) = point[2]^(-2.0)
scg33i(point::Vector{<: Number}) = point[2]^(-2.0) * csc(point[3])

scmetric = Metric([scg00, scg11, scg22, scg33], [scg00i, scg11i, scg22i, scg33i])

t = 0.0
r = 6.0M
θ = π/2
ϕ = 0.0

point  = Vector([t, r, θ, ϕ]);
u_form = Vector([1.0, 0.0, 0.0, 0.0])

# Tests
@test typeof(point) <: AbstractVector
@test eltype(point) <: Number

@test dim(scmetric) == 4

gμν, igμν = evaluate(scmetric, point)

@test igμν*gμν ≈ diagm(ones(dim(scmetric)))

scfactor = (1.0 - (2.0M/r))

testH = (scfactor)^(1//2) * (scfactor*u_form[2]^2 + (u_form[3] / r)^2 + (u_form[4] / (r * sin(θ)))^2 + 1.0)^(1//2)

@test hamiltonian(gμν, point, u_form) ≈ testH
