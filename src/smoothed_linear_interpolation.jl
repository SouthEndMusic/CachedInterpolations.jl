"""
    SmoothedLinearInterpolation(u, t; λ = 0.25, extrapolate = false)

The method of interpolating between the data points using a linear polynomial, with an anterval around the
corner points being replaced by a smooth spline section.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `λ`: The relative size of the spline interval. The interval extents a fraction `λ/2` towards
    the neighbouring time points.
"""
struct SmoothedLinearInterpolation{uType, tType, λType <: Real, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    cache::SmoothedLinearInterpolationCache{uType, tType, λType}
    extrapolate::Bool
    function SmoothedLinearInterpolation(u, t, cache, λ, extrapolate)
        return new{typeof(u), typeof(t), typeof(λ), eltype(u)}(u, t, cache, extrapolate)
    end
end

function SmoothedLinearInterpolation(
    u,
    t;
    λ = 0.25,
    extrapolate::Bool = false,
)::SmoothedLinearInterpolation
    u, t = munge_data(u, t)
    # Make sure the parameter λ is in the right range
    @assert 0 <= λ <= 1 "The parameter λ must be in the interval [0,1], got $λ."
    cache = SmoothedLinearInterpolationCache(u, t, λ)
    return SmoothedLinearInterpolation(u, t, cache, λ, extrapolate)
end

"""
    S(A, t, idx)

Compute the spline parameter `s` from from the time `t`.

    ## Arguments

    - `A`: The `SmoothedLinearInterpolation` object
    - `t`: The time point
    - `idx`: The index indicating which spline section
"""
function S(A::SmoothedLinearInterpolation, t, idx)
    (; Δt, ΔΔt, t_tilde, λ) = A.cache
    Δtᵢ = Δt[idx]
    ΔΔtᵢ = ΔΔt[idx]
    tdiff = t - t_tilde[2 * idx - 1]
    @assert tdiff >= 0

    if isapprox(ΔΔtᵢ, 0; atol = 1e-5)
        # Degenerate case Δtᵢ₊₁ ≈ Δtᵢ
        1 / λ * tdiff / Δtᵢ
    else
        (-Δtᵢ + sqrt(Δtᵢ^2 + 2 * ΔΔtᵢ * tdiff / λ)) / ΔΔtᵢ
    end
end

"""
    S(A, t, idx)

Compute the spline value `u` at the time `t`.

    ## Arguments

    - `A`: The `SmoothedLinearInterpolation` object
    - `t`: The time point
    - `idx`: The index indicating which spline section
"""
function U(A::SmoothedLinearInterpolation, t, idx)
    s = S(A, t, idx)
    return U_s(A, s, idx)
end

"""
    S(A, t, idx)

Compute the spline value `u` from the spline parameter `s`.

    ## Arguments

    - `A`: The `SmoothedLinearInterpolation` object
    - `s`: The spline parameter value
    - `idx`: The index indicating which spline section
"""
function U_s(A::AbstractInterpolation, s, idx)
    (; Δu, ΔΔu, u_tilde, λ) = A.cache
    Δuᵢ = Δu[idx]
    ΔΔuᵢ = ΔΔu[idx]
    return λ / 2 * ΔΔuᵢ * s^2 + λ * Δuᵢ * s + u_tilde[2 * idx - 1]
end

function DataInterpolations._interpolate(
    A::SmoothedLinearInterpolation{<:AbstractVector},
    t::Number,
    iguess,
)
    (; u, t_tilde, linear_slope) = A.cache

    # idx of smallest idx such that A.t[idx] >= t
    idx = searchsortedfirstcorrelated(A.t, t, iguess)

    if idx == 1
        # Linear extrapolation for t < A.t[1]
        u[1] + linear_slope[1] * (t - A.t[1])
    elseif idx == length(u) + 1
        # Linear extrapolation for t > A.t[end]
        u[end] + linear_slope[end] * (t - A.t[end])
    else
        # Interpolation
        if t < t_tilde[2 * idx - 2]
            U(A, t, idx - 1)
        elseif t > t_tilde[2 * idx - 1]
            U(A, t, idx)
        else
            # Linear interpolation
            u[idx] + linear_slope[idx] * (t - A.t[idx])
        end
    end
end
