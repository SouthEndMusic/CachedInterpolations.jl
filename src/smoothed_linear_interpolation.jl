struct CLinearInterpolation{uType, tType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    cache::LinearInterpolationCache{uType}
    extrapolate::Bool
    function CLinearInterpolation(u, t, cache, extrapolate)
        new{typeof(u), typeof(t), eltype(u)}(u, t, cache, extrapolate)
    end
end

function CLinearInterpolation(u, t; extrapolate = false)::CLinearInterpolation
    u, t = munge_data(u, t)
    cache = LinearInterpolationCache(u, t)
    return CLinearInterpolation(u, t, cache, extrapolate)
end

function DataInterpolations._interpolate(
    A::CLinearInterpolation{<:AbstractVector},
    t::Number,
    iguess,
)
    (; slope, idx_prev) = A.cache

    # idx of smallest idx such that A.t[idx] >= t
    idx = searchsortedfirstcorrelated(A.t, t, idx_prev[])
    idx = min(idx, length(slope))
    idx_prev[] = idx

    return A.u[idx] + slope[idx] * (t - A.t[idx])
end

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
    linear_itp::LinearInterpolation{uType, tType, T}
    extrapolate::Bool
    function SmoothedLinearInterpolation(u, t, cache, λ, linear_itp, extrapolate)
        return new{typeof(u), typeof(t), typeof(λ), eltype(u)}(
            u,
            t,
            cache,
            linear_itp,
            extrapolate,
        )
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
    linear_itp = LinearInterpolation(u, t; extrapolate)
    return SmoothedLinearInterpolation(u, t, cache, λ, linear_itp, extrapolate)
end

function DataInterpolations._interpolate(
    A::SmoothedLinearInterpolation{<:AbstractVector},
    t::Number,
    iguess,
)
    (; u, t_tilde, idx_prev) = A.cache

    # idx of smallest idx such that A.t[idx] >= t
    idx = searchsortedfirstcorrelated(A.t, t, idx_prev[])
    idx_prev[] = idx

    if idx == 1 || idx == length(u) + 1
        # Linear extrapolation
        DataInterpolations._interpolate(A.linear_itp, t, idx)[1]
    else
        # Interpolation
        if t < t_tilde[2 * idx - 2]
            U(A, t, idx - 1)
        elseif t > t_tilde[2 * idx - 1]
            U(A, t, idx)
        else
            # Linear interpolation
            DataInterpolations._interpolate(A.linear_itp, t, idx)[1]
        end
    end
end
