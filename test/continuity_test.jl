using SmoothInterpolation
using Random
using ForwardDiff

@testset "CSmoothedLinearInterpolation degenerate" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)

    u₋ = @. itp(itp.cache.t_tilde - ε)
    u₊ = @. itp(itp.cache.t_tilde + ε)

    @test u₋ ≈ u₊ atol = 1e-4
end

@testset "CSmoothedLinearInterpolation non-degenerate" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = cumsum(rand(5) .+ (1:5))
    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)

    u₋ = @. itp(itp.cache.t_tilde - ε)
    u₊ = @. itp(itp.cache.t_tilde + ε)

    @test u₋ ≈ u₊ atol = 1e-4
end

@testset "CSmoothedLinearInterpolation degenerate derivative" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)

    du₋ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .- ε)
    du₊ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .+ ε)

    @test du₋ ≈ du₊ atol = 1e-4
end

@testset "CSmoothedLinearInterpolation non-degenerate derivative" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = cumsum(rand(5) .+ (1:5))
    itp = CSmoothedLinearInterpolation(u, t; extrapolate = true)

    du₋ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .- ε)
    du₊ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .+ ε)

    @test du₋ ≈ du₊ atol = 1e-4
end