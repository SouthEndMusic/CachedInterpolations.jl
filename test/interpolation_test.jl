using DataInterpolations
using CachedInterpolations
using Random

@testset "CLinearInterpolation" begin
    Random.seed!(42)

    u = cumsum(rand(5))
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = LinearInterpolation(u, t)
    citp = CLinearInterpolation(u, t)

    t_eval = 1.5:0.3:5.0
    @test citp.(t_eval) ≈ itp.(t_eval)
end

@testset "CSmoothedLinearInterpolation degenerate" begin
    Random.seed!(1)

    u = cumsum(rand(5))
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = CSmoothedLinearInterpolation(u, t)

    @test all(itp.cache.ΔΔt .≈ 0)
    @test itp.(1.5:0.3:5.0) ≈ [
        0.24799,
        0.35276,
        0.49293,
        0.70214,
        0.91179,
        1.11923,
        1.30991,
        1.49839,
        1.68723,
        1.93269,
        2.20716,
        2.48164,
    ] atol = 1e-4
    @test_nowarn string(itp.cache)
end

@testset "CSmoothedLinearInterpolation non-degenerate" begin
    Random.seed!(2)
    u = cumsum(rand(5))
    t = cumsum(rand(5) .+ (1:5))
    itp = CSmoothedLinearInterpolation(u, t)

    @test !any(itp.cache.ΔΔt[2:(end - 1)] .≈ 0)
    @test itp.(t[1]:1.2:t[end]) ≈ [
        0.00225,
        0.30983,
        0.61718,
        0.88465,
        1.14094,
        1.39383,
        1.60492,
        1.80825,
        2.01157,
        2.19797,
        2.27307,
        2.32361,
        2.37415,
        2.42469,
    ] atol = 1e-4
    @test_nowarn string(itp)
end

@testset "CLinearInterpolationIntInv" begin
    Random.seed!(9)
    u = rand(5)
    # Add degenerate case of constant u
    push!(u, u[end])
    t = cumsum(rand(6))
    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)
    itp = CLinearInterpolation(itp)
    itp_int_inv = invert_integral(itp)
    t_eval = range(t[1], t[end]; length = 200)
    u_int_eval = DataInterpolations.integral.(Ref(itp), t_eval)
    @test t_eval ≈ itp_int_inv.(u_int_eval)
    @test_nowarn string(itp_int_inv.cache)
end

@testset "CSmoothedLinearInterpolationIntInv" begin
    Random.seed!(9)
    u = rand(5)
    # Add degenerate case of constant u
    push!(u, u[end])
    t = cumsum(rand(6))
    itp = CSmoothedLinearInterpolation(u, t)
    itp_int_inv = invert_integral(itp)
    t_eval = range(t[1], t[end]; length = 200)
    u_int_eval = DataInterpolations.integral.(Ref(itp), t_eval)
    @test t_eval ≈ itp_int_inv.(u_int_eval)
    @test_nowarn string(itp_int_inv.cache_integration)
end
