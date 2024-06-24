abstract type AbstractCache{T} end

struct CLinearInterpolationCache{uType, T} <: AbstractCache{T}
    slope::uType
    idx_prev::Base.RefValue{Int}
    function CLinearInterpolationCache(slope, idx_prev::Base.RefValue{Int})
        new{typeof(slope), eltype(slope)}(slope, idx_prev)
    end
end

function CLinearInterpolationCache(u, t)
    Δu = diff(u)
    Δt = diff(t)
    slope = Δu ./ Δt
    pushfirst!(slope, first(slope))
    return CLinearInterpolationCache(slope, Ref(1))
end

"""
The cache object for CLinearInterpolationIntInv computations.
"""
struct CLinearInterpolationIntInvCache{uType, T} <: AbstractCache{T}
    u::uType
    slope::uType
    degenerate_slope::Vector{Bool}
    idx_prev::Base.RefValue{Int}
end

function CLinearInterpolationIntInvCache(u, t)
    Δu = diff(u)
    Δt = diff(t)
    slope = Δu ./ Δt
    degenerate_slope = collect(isapprox.(slope, 0, atol = 1e-5))
    return CLinearInterpolationIntInvCache{typeof(u), eltype(u)}(
        u,
        slope,
        degenerate_slope,
        Ref(1),
    )
end

"""
The cache object for CSmoothedLinearInterpolation computations.
"""
struct CSmoothedLinearInterpolationCache{uType, tType, λType <: Number, T} <:
       AbstractCache{T}
    u::uType
    t::tType
    Δu::uType
    Δt::tType
    ΔΔu::uType
    ΔΔt::tType
    u_tilde::uType
    t_tilde::tType
    linear_slope::uType
    # Whether ΔΔt is sufficiently close to 0
    degenerate_ΔΔt::Vector{Bool}
    λ::λType
    idx_prev::Base.RefValue{Int}
end

function CSmoothedLinearInterpolationCache(u, t, λ)::CSmoothedLinearInterpolationCache
    Δu = diff(u)
    Δt = diff(t)
    @assert !any(iszero.(Δt))
    pushfirst!(Δt, first(Δt))
    push!(Δt, last(Δt))
    pushfirst!(Δu, first(Δu))
    push!(Δu, last(Δu))
    ΔΔu = diff(Δu)
    ΔΔt = diff(Δt)
    u_tilde = get_spline_ends(u, Δu, λ)
    t_tilde = get_spline_ends(t, Δt, λ)
    linear_slope = Δu ./ Δt
    # Whether ΔΔt is sufficiently close to 0
    degenerate_ΔΔt = collect(isapprox.(ΔΔt, 0, atol = 1e-5))
    return CSmoothedLinearInterpolationCache{typeof(u), typeof(t), typeof(λ), eltype(u)}(
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
        Ref(1),
    )
end

"""
The cache object for CSmoothedLinearInterpolationIntInv computations.
"""
struct CSmoothedLinearInterpolationIntInvCache{uType, T} <: AbstractCache{T}
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
    # Whether Δu is sufficiently close to 0
    degenerate_Δu::Vector{Bool}
    idx_prev::Base.RefValue{Int}
end

function CSmoothedLinearInterpolationIntInvCache(A)
    coeffs =
        hcat([collect(get_quartic_coefficients(A, idx)) for idx in eachindex(A.cache.t)]...)
    c4, c3, c2, c1 = collect.(eachrow(coeffs))
    # The degree is 5 minus the index of the first (≈) nonzero coefficient
    degree = 5 .- findfirst.(coef -> !isapprox(coef, 0; atol = 1e-5), eachcol(coeffs))

    p = p_coeff.(c4, c3, c2)
    q = q_coeff.(c4, c3, c2, c1)

    # Whether Δu is sufficiently close to 0
    degenerate_Δu = collect(isapprox.(A.cache.Δu, 0, atol = 1e-5))

    return CSmoothedLinearInterpolationIntInvCache{typeof(A.u), eltype(A.u)}(
        degree,
        c4,
        c3,
        c2,
        c1,
        p,
        q,
        degenerate_Δu,
        Ref(1),
    )
end
