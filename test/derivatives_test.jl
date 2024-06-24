using Random
using DataInterpolations
using CachedInterpolations
using ForwardDiff

@testset "CSmoothedLinearInterpolation" begin
    Random.seed!(10)

    t = cumsum(rand(10))
    u = rand(10)
    t_eval = (t[1] - 1):0.01:(t[end] + 1)

    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)
    u_deriv_eval = DataInterpolations.derivative.(Ref(itp), t_eval)
    u_deriv_check = ForwardDiff.derivative.(Ref(itp), t_eval)

    @test u_deriv_eval ≈ u_deriv_check
end

@testset "CLinearInterpolationIntInv" begin
    Random.seed!(10)

    t = cumsum(rand(10))
    u = rand(10)

    itp = LinearInterpolation(u, t; extrapolate = true)
    itp_int_inv = invert_integral(itp)
    u_int_eval = itp_int_inv.t[1]:0.01:(itp_int_inv.t[end] + 1)

    t_deriv_eval = DataInterpolations.derivative.(Ref(itp_int_inv), u_int_eval)
    t_deriv_check = ForwardDiff.derivative.(Ref(itp_int_inv), u_int_eval)

    @test t_deriv_eval ≈ t_deriv_check
end

@testset "CSmoothedLinearInterpolationIntInv" begin
    Random.seed!(10)

    t = cumsum(rand(10))
    u = rand(10)

    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)
    itp_int_inv = invert_integral(itp)
    u_int_eval = itp_int_inv.t[1]:0.01:(itp_int_inv.t[end] + 1)

    t_deriv_eval = DataInterpolations.derivative.(Ref(itp_int_inv), u_int_eval)
    t_deriv_check = ForwardDiff.derivative.(Ref(itp_int_inv), u_int_eval)

    @test t_deriv_eval ≈ t_deriv_check
end
