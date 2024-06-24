using CachedInterpolations: iterate_roots, p_coeff, q_coeff

@testset "Polynomial solving" begin
    p = 0.0
    q = 0.0

    # Degree 1
    c4 = 0.0
    c3 = 0.0
    c2 = 0.0
    c1 = 2.0
    c0 = 4.0
    roots = collect(iterate_roots(1, c4, c3, c2, c1, c0, p, q))
    @test length(roots) == 1
    @test only(roots) ≈ -2.0

    # Degree 2
    c4 = 0.0
    c3 = 0.0
    c2 = 1.0
    c1 = -1.0
    c0 = -1.0
    roots = collect(iterate_roots(2, c4, c3, c2, c1, c0, p, q))
    @test length(roots) == 2
    ϕ = (sqrt(5) + 1) / 2
    @test roots ≈ [1 - ϕ, ϕ]

    # Degree 3
    c4 = 0.0
    c3 = 3.0
    c2 = -63.0
    c1 = 429.0
    c0 = -945.0
    roots = collect(iterate_roots(3, c4, c3, c2, c1, c0, p, q))
    @test roots ≈ Float64[5, 9, 7]

    # Degree 4
    c4 = 1.0
    c3 = -10.0
    c2 = 35.0
    c1 = -50.0
    c0 = 24.0
    p = p_coeff(c4, c3, c2)
    q = q_coeff(c4, c3, c2, c1)
    roots = collect(iterate_roots(4, c4, c3, c2, c1, c0, p, q))
    @test roots ≈ Float64[2, 4, 3, 1]
end