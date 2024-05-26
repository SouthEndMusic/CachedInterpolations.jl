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

function SmoothedLinearInterpolationIntInv(
    A::SmoothedLinearInterpolation,
)::SmoothedLinearInterpolationIntInv
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
    (; a, b, c, d, p, q) = cache_integration

    # idx of smallest idx such that A.t[idx] >= t
    # Note that A.t denotes integrated values
    idx = searchsortedfirstcorrelated(A.t, V, iguess)
    @show V

    if idx == 1
        @assert V >= 0 "Cannot invert intagral for negative input."
        idx = 2
    end

    if idx == 2
        # First half spline section
        # which is linear
        Vdiff = (V - t[1])
        @assert Vdiff >= 0
        u[1] +
        (-cache.u[1] + sqrt(cache.u[1]^2 + 2 * cache.linear_slope[1] * Vdiff)) /
        cache.linear_slope[1]
    elseif idx == n_points + 1
        # Extrapolation
        Vdiff = (V - t[end])
        @assert Vdiff >= 0
        u[end] +
        (-cache.u[end] + sqrt(cache.u[end]^2 + 2 * cache.linear_slope[end] * Vdiff)) /
        cache.linear_slope[end]
    elseif idx % 2 == 0
        Vdiff = (V - A.t[idx - 1])
        @assert Vdiff >= 0

        i = Int(idx // 2)
        aᵢ = a[i]
        bᵢ = b[i]
        cᵢ = c[i]
        dᵢ = d[i]
        pᵢ = p[i]
        qᵢ = q[i]

        Δ₀ = Complex(cᵢ^2 - 3 * bᵢ * dᵢ - 12 * aᵢ * Vdiff)
        Δ₁ = Complex(
            2 * cᵢ^3 - 9 * bᵢ * cᵢ * dᵢ - 27 * bᵢ^2 * Vdiff +
            27 * aᵢ * dᵢ^2 +
            72 * aᵢ * cᵢ * Vdiff,
        )
        Q = ((Δ₁ + sqrt(Δ₁^2 - 4 * Δ₀^3)) / 2)^(1 / 3)
        S = sqrt(-2 * pᵢ / 3 + (Q + Δ₀ / Q) / (3 * aᵢ)) / 2

        root1 = sqrt(-4 * S^2 - 2 * pᵢ + qᵢ / S)
        root2 = sqrt(-4 * S^2 - 2 * pᵢ - qᵢ / S)

        # Check the 4 possible roots for being valid;
        # real and in [0,1]
        s1 = -bᵢ / (4 * aᵢ) - S + root1 / 2
        if valid(s1)
            return U_s(A, real(s1), i)
        end

        s2 = -bᵢ / (4 * aᵢ) - S - root1 / 2
        if valid(s2)
            return U_s(A, real(s2), i)
        end

        s3 = -bᵢ / (4 * aᵢ) + S + root2 / 2
        if valid(s3)
            return U_s(A, real(s3), i)
        end

        s4 = -bᵢ / (4 * aᵢ) + S - root2 / 2
        if valid(s4)
            return U_s(A, real(s4), i)
        end

        @show aᵢ
        error("No valid root found, got $([s1,s2,s3,s4]).")
    else
        # Linear section of SmoothedLinearInterpolation
        Vdiff = (V - A.t[idx - 1])
        @assert Vdiff >= 0

        # TODO: Consider degenerate case
        # where Δu[i+1] = 0
        i = Int((idx - 1) // 2)
        Δuᵢ₊₁ = cache.Δu[i + 1]
        Δtᵢ₊₁ = cache.Δt[i + 1]
        u_frac = cache.u[i] / Δuᵢ₊₁
        λ = cache.λ
        root = sqrt(u_frac^2 + λ * (u_frac + λ / 4) + 2 * Vdiff / (Δtᵢ₊₁ * Δuᵢ₊₁))
        cache.t[i] + (-u_frac + sign(u_frac) * root) * Δtᵢ₊₁
    end
end
