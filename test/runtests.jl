using GeodesicIntegration
using Base.Test

@time @testset "Vector and Dual Vector tests" begin include("vector_test.jl") end
@time @testset "Numerical algorithm tests" begin include("math_test.jl") end
@time @testset "Schwarzschild metric tests" begin include("schwarzschild_test.jl") end
