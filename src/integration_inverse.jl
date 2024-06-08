"""
    LinearInterpolationIntInv(A::SmoothedLinearInterpolation)

Inverting the integral of a LinearInterpolation object if possible. The `A.u` must be non-negative.

## Arguments

  - A The LinearInterpolation object whose integral is inverted.
"""
struct LinearInterpolationIntInv{uType, tType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    cache::LinearInterpolationIntInvCache{uType}
    extrapolate::Bool
    function LinearInterpolationIntInv(u, t, cache, extrapolate)
        return new{typeof(u), typeof(t), eltype(u)}(u, t, cache, extrapolate)
    end
end

"""
    Invert the integral of a LinearInterpolation object, which yields 
    LinearInterpolationIntInv object.
"""
function invert_integral(A::LinearInterpolation)::LinearInterpolationIntInv
    @assert all(A.u .>= 0) "Inverting the integral is only supported for non-negative LinearInterpolation."
    t = DataInterpolations.integral.(Ref(A), A.t)
    cache = LinearInterpolationIntInvCache(A.u, A.t)
    return LinearInterpolationIntInv(A.t, t, cache, A.extrapolate)
end

function DataInterpolations._interpolate(
    A::LinearInterpolationIntInv{<:AbstractVector},
    V::Number,
    iguess,
)
    (; cache) = A

    # idx of smallest idx such that A.t[idx] >= V
    # Note that A.t denotes integrated values
    idx = searchsortedfirstcorrelated(A.t, V, iguess)

    if idx == 1
        @assert V >= 0 "Cannot invert integral for negative input."
        idx = 2
    end

    Vdiff = (V - A.t[idx - 1])
    @assert Vdiff >= 0 "Vdiff must be non_negative, got V = $V, Vdiff = $Vdiff, idx = $idx"

    t_prev = A.u[idx - 1]
    idx = min(idx, length(A.u))

    i = idx - 1

    if cache.degenerate_Δu[i]
        # Special case when LinearInterpolation is (near) constant
        t_prev + Vdiff / cache.u[idx]
    else
        t_prev +
        (-cache.u[i] + sqrt(cache.u[i]^2 + 2 * cache.slope[i] * Vdiff)) / cache.slope[i]
    end
end

"""
    SmoothedLinearInterpolationIntInv(A::SmoothedLinearInterpolation)

Inverting the integral of a SmoothedLinearInterpolation object if possible. The `A.u` must be non-negative.

## Arguments

  - A The SmoothedLinearInterpolation object whose integral is inverted.
"""
struct SmoothedLinearInterpolationIntInv{uType, tType, λType <: Real, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    cache::SmoothedLinearInterpolationCache{uType, tType, λType}
    cache_integration::SmoothedLinearInterpolationIntInvCache{uType}
    extrapolate::Bool
    function SmoothedLinearInterpolationIntInv(u, t, cache, cache_int, λ, extrapolate)
        return new{typeof(u), typeof(t), typeof(λ), eltype(u)}(
            u,
            t,
            cache,
            cache_int,
            extrapolate,
        )
    end
end

"""
    Invert the integral of a SmoothedLinearInterpolation object, which yields 
    SmoothedLinearInterpolationIntInv object.
"""
function invert_integral(A::SmoothedLinearInterpolation)::SmoothedLinearInterpolationIntInv
    @assert all(A.u .>= 0) "Inverting the integral is only supported for non-negative SmoothedLinearInterpolation."
    (; cache, extrapolate) = A
    t = DataInterpolations.integral.(Ref(A), cache.t_tilde)
    u = cache.t_tilde
    cache_int = SmoothedLinearInterpolationIntInvCache(A)
    return SmoothedLinearInterpolationIntInv(u, t, cache, cache_int, cache.λ, extrapolate)
end

function DataInterpolations._interpolate(
    A::SmoothedLinearInterpolationIntInv{<:AbstractVector},
    V::Number,
    iguess,
)
    n_points = length(A.t)
    (; u, t, cache, cache_integration) = A
    (; degree, c4, c3, c2, c1, p, q, degenerate_Δu) = cache_integration

    # idx of smallest idx such that A.t[idx] >= V
    # Note that A.t denotes integrated values
    idx = searchsortedfirstcorrelated(A.t, V, iguess)

    if idx == 1
        @assert V >= 0 "Cannot invert integral for negative input."
        idx = 2
    end

    if idx == 2
        # First half spline section
        # which is linear
        Vdiff = (V - t[1])
        @assert Vdiff >= 0
        if isapprox(cache.linear_slope[1], 0; atol = 1e-5)
            u[1] + Vdiff / cache.u[1]
        else
            u[1] +
            (-cache.u[1] + sqrt(cache.u[1]^2 + 2 * cache.linear_slope[1] * Vdiff)) /
            cache.linear_slope[1]
        end
    elseif idx == n_points + 1
        # Extrapolation
        Vdiff = (V - t[end])
        @assert Vdiff >= 0
        if isapprox(cache.linear_slope[end], 0; atol = 1e-5)
            u[end] + Vdiff / cache.u[end]
        else
            u[end] +
            (-cache.u[end] + sqrt(cache.u[end]^2 + 2 * cache.linear_slope[end] * Vdiff)) /
            cache.linear_slope[end]
        end
    elseif idx % 2 == 0
        Vdiff = (V - A.t[idx - 1])
        @assert Vdiff >= 0

        i = idx ÷ 2
        c4ᵢ = c4[i]
        c3ᵢ = c3[i]
        c2ᵢ = c2[i]
        c1ᵢ = c1[i]
        c0 = -Vdiff
        pᵢ = p[i]
        qᵢ = q[i]
        degᵢ = degree[i]

        # Check the 4 possible roots for being valid;
        # real and in [0,1]
        root_iterator = iterate_roots(degᵢ, c4ᵢ, c3ᵢ, c2ᵢ, c1ᵢ, c0, pᵢ, qᵢ)
        for s in root_iterator
            if valid(s)
                return T_s(A, real(s), i)
            end
        end

        error("No valid root found, got $(collect(root_iterator)) for V = $V.")
    else
        # Linear section of SmoothedLinearInterpolation
        Vdiff = (V - A.t[idx - 1])
        @assert Vdiff >= 0
        i = (idx - 1) ÷ 2

        if degenerate_Δu[i + 1]
            # Special case when SmoothedLinearInterpolation is (near) constant
            A.u[idx - 1] + Vdiff / cache.u[i]
        else
            Δuᵢ₊₁ = cache.Δu[i + 1]
            Δtᵢ₊₁ = cache.Δt[i + 1]
            u_frac = cache.u[i] / Δuᵢ₊₁
            λ = cache.λ
            root = sqrt(u_frac^2 + λ * (u_frac + λ / 4) + 2 * Vdiff / (Δtᵢ₊₁ * Δuᵢ₊₁))
            cache.t[i] + (-u_frac + sign(u_frac) * root) * Δtᵢ₊₁
        end
    end
end