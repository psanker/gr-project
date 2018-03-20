using GeodesicIntegration
using Base.Test

@time @testset "Numerical algorithm tests" begin include("math_test.jl") end
@time @testset "Schwarzschild metric tests" begin include("schwarzschild_test.jl") end
