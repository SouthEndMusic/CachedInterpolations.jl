using SmoothInterpolation
using Random

@testset "SmoothedLinearInterpolation degenerate" begin
    Random.seed!(1)

    u = cumsum(rand(5))
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = SmoothedLinearInterpolation(u, t)

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
end

@testset "SmoothedLinearInterpolation non-degenerate" begin
    Random.seed!(2)
    u = cumsum(rand(5))
    t = cumsum(rand(5) .+ (1:5))
    itp = SmoothedLinearInterpolation(u, t)

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
end