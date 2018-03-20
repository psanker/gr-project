using GeodesicIntegration, Base.Test

point = 0.0

@test ∂(exp, point) ≈ 1.0 # Direct numerical computation

pts = randn(100)
res = zeros(length(pts))

for ind ∈ eachindex(pts)
    res[ind] = ∂(exp, pts[ind])
end

@test res ≈ exp.(pts)

# Verify gradients on specific variables
α = 2.0
β = 3.0
f(x::Vector{<: Number}) = exp(β * x[1]) * exp(α * x[2])

point = [4.5, 3.0, 29., 923.]

@test ∂(f, point) ≈ [β * f(point), α * f(point), 0.0, 0.0]
@test ∂(f, point; set=[2, 3, 4]) ≈ [0.0, α * f(point), 0.0, 0.0]
