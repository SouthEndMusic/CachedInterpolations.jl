"""
The cache object for SmoothedLinearInterpolation computations.
"""
struct SmoothedLinearInterpolationCache{uType, tType, λType <: Number}
    u::uType
    t::tType
    Δu::uType
    Δt::tType
    ΔΔu::uType
    ΔΔt::tType
    u_tilde::uType
    t_tilde::tType
    linear_slope::uType
    degenerate_ΔΔt::BitVector
    λ::λType
end

function SmoothedLinearInterpolationCache(u, t, λ)::SmoothedLinearInterpolationCache
    Δu = diff(u)
    Δt = diff(t)
    @assert !any(iszero.(Δt))
    pushfirst!(Δt, Δt[1])
    push!(Δt, Δt[end])
    pushfirst!(Δu, Δu[1])
    push!(Δu, Δu[end])
    ΔΔu = diff(Δu)
    ΔΔt = diff(Δt)
    u_tilde = get_spline_ends(u, Δu, λ)
    t_tilde = get_spline_ends(t, Δt, λ)
    linear_slope = Δu ./ Δt
    degenerate_ΔΔt = BitVector(isapprox.(ΔΔt, 0, atol = 1e-5))
    return SmoothedLinearInterpolationCache(
        u,
        t,
        Δu,
        Δt,
        ΔΔu,
        ΔΔt,
        u_tilde,
        t_tilde,
        linear_slope,
        degenerate_ΔΔt,
        λ,
    )
end

"""
The cache object for SmoothedLinearInterpolationIntInv computations.
"""
struct SmoothedLinearInterpolationIntInvCache{uType}
    # The degree of the polynomial whose roots need to be found
    degree::Vector{Int}
    # Quartic polynomial coefficients
    c4::uType
    c3::uType
    c2::uType
    c1::uType
    # Coefficients of depressed quartic
    p::uType
    q::uType
    # Degeneracy case
    degenerate_Δu::BitVector
end

function SmoothedLinearInterpolationIntInvCache(A)
    coeffs =
        hcat([collect(get_quartic_coefficients(A, idx)) for idx in eachindex(A.cache.t)]...)
    c4, c3, c2, c1 = collect.(eachrow(coeffs))
    # The degree is 5 minus the index of the first (≈) nonzero coefficient
    degree = 5 .- findfirst.(coef -> !isapprox(coef, 0; atol = 1e-5), eachcol(coeffs))

    p = p_coeff.(c4, c3, c2)
    q = q_coeff.(c4, c3, c2, c1)

    degenerate_Δu = isapprox.(A.cache.Δu, 0, atol = 1e-5)

    return SmoothedLinearInterpolationIntInvCache(
        degree,
        c4,
        c3,
        c2,
        c1,
        p,
        q,
        degenerate_Δu,
    )
end
