using CachedInterpolations

@testset "CSmoothedLinearInterpolation" begin
    u = [1.0, 2.0, 3.0, 4.0]
    t = [1.0, 2.0, 3.0, 4.0]

    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)
    @test itp(0.0) ≈ 0.0
    @test itp(5.0) ≈ 5.0

    u = zeros(5)
    t = [1.0, 2.0, 3.0, 4.0, 5.0]

    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)
    @test itp(0.0) ≈ 0.0
    @test itp(5.0) ≈ 0.0
end