using SmoothInterpolation
using Random
using ForwardDiff

@testset "SmoothedLinearInterpolation degenerate" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)

    u₋ = @. itp(itp.cache.t_tilde - ε)
    u₊ = @. itp(itp.cache.t_tilde + ε)

    @test u₋ ≈ u₊ atol = 1e-4
end

@testset "SmoothedLinearInterpolation non-degenerate" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = cumsum(rand(5) .+ (1:5))
    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)

    u₋ = @. itp(itp.cache.t_tilde - ε)
    u₊ = @. itp(itp.cache.t_tilde + ε)

    @test u₋ ≈ u₊ atol = 1e-4
end

@testset "SmoothedLinearInterpolation degenerate derivative" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)

    du₋ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .- ε)
    du₊ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .+ ε)

    @test du₋ ≈ du₊ atol = 1e-4
end

@testset "SmoothedLinearInterpolation non-degenerate derivative" begin
    Random.seed!(1)
    ε = 1e-5

    u = cumsum(rand(5))
    t = cumsum(rand(5) .+ (1:5))
    itp = SmoothedLinearInterpolation(u, t; extrapolate = true)

    du₋ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .- ε)
    du₊ = ForwardDiff.derivative.(Ref(itp), itp.cache.t_tilde .+ ε)

    @test du₋ ≈ du₊ atol = 1e-4
end