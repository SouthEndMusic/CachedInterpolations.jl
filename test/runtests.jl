using Test, InterpolationConcepts
using Random: seed!

seed!(1)
u = rand(10)
t = cumsum(rand(10))
t_eval = cumsum(rand(15))

@testset "SomeStuff" begin
    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)
    itp.(t_eval)
end
