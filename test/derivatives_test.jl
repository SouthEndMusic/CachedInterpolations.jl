using Random
using DataInterpolations
using SmoothInterpolation
using ForwardDiff

@testset "SmoothedLinearInterpolation" begin
    Random.seed!(10)

    t = cumsum(rand(10))
    u = rand(10)
    t_eval = (t[1] - 1):0.01:(t[end] + 1)

    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)
    u_deriv_eval = DataInterpolations.derivative.(Ref(itp), t_eval)
    u_deriv_check = ForwardDiff.derivative.(Ref(itp), t_eval)

    @test u_deriv_eval ≈ u_deriv_check
end