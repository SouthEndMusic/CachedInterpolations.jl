using DataInterpolations
using SmoothInterpolation
using Random
using ForwardDiff

@testset "SmoothedLinearInterpolation integration outcome" begin
    Random.seed!(1)

    u = rand(5)
    t = cumsum(rand(5))
    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)

    # With extrapolation
    t_eval = t[1]:0.1:(t[end] + 1)
    u_int = DataInterpolations.integral.(Ref(itp), t_eval)

    # Numerical integration
    u_eval = itp.(t_eval)
    u_int_num = 0.5 * 0.1 * (u_eval[2:end] + u_eval[1:(end - 1)])
    u_int_num = cumsum(u_int_num)
    pushfirst!(u_int_num, 0.0)

    @test u_int ≈ u_int_num rtol = 1e-3
end

@testset "SmoothedLinearInterpolation integration derivative" begin
    Random.seed!(1)
    t = cumsum(rand(5))

    u = rand(5)
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)

    # With extrapolation
    t_eval = t[1]:0.1:(t[end] + 1)
    u_eval = itp.(t_eval)

    integral = t -> DataInterpolations.integral(itp, t)
    u_eval_ = ForwardDiff.derivative.(Ref(integral), t_eval)

    # Automatic derivative at t = itp.t[1] fails (#25)
    @test u_eval[2:end] ≈ u_eval_[2:end]
end