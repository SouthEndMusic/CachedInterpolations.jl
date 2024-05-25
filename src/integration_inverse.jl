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
    t = integrate_sections(cache)
    u = cache.t_tilde
    cache_int = SmoothedLinearInterpolationIntInvCache(A)
    return SmoothedLinearInterpolationIntInv(u, t, cache, cache_int, cache.λ, extrapolate)
end

function integrate_sections(cache::SmoothedLinearInterpolationCache)
    (; u, Δu, Δt, u_tilde, t_tilde, λ) = cache
    n_points = length(t_tilde)
    U = zeros(n_points)
    for j in eachindex(t_tilde)
        if j == 1
            continue
        elseif j == 2
            U[j] = 0.5 * (u_tilde[2] + u_tilde[1]) * (t_tilde[2] - t_tilde[1])
        elseif j == n_points
            U[j] =
                0.5 * (u_tilde[end] + u_tilde[end - 1]) * (t_tilde[end] - t_tilde[end - 1])
        elseif j % 2 == 0
            i = Int(j // 2)
            Δuᵢ = Δu[i]
            Δuᵢ₊₁ = Δu[i + 1]
            Δtᵢ = Δt[i]
            Δtᵢ₊₁ = Δt[i + 1]
            U[j] =
                λ^2 / 24 * (Δtᵢ * (-3 * Δuᵢ + Δuᵢ₊₁) + Δtᵢ₊₁ * (-Δuᵢ + 3 * Δuᵢ₊₁)) +
                λ / 2 * (Δtᵢ + Δtᵢ₊₁) * u[i]
        else
            U[j] = 0.5 * (u_tilde[j] + u_tilde[j - 1]) * (t_tilde[j] - t_tilde[j - 1])
        end
    end
    U = cumsum(U)
    return U
end

function DataInterpolations._interpolate(
    A::SmoothedLinearInterpolationIntInv{<:AbstractVector},
    V::Number,
    iguess,
)
    n_points = length(A.t)
    (; u, t, Δu, Δt, λ) = A.cache
    (; a, b, c, d, p, q) = A.cache_integration

    # idx of smallest idx such that A.t[idx] >= t
    # Note that A.t denotes integrated values
    idx = searchsortedfirstcorrelated(A.t, V, iguess)

    if idx == 1
        0
    elseif idx == 2
        0
    elseif idx == n_points + 1
        0
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
            2 * cᵢ^3 +
            9 * bᵢ * cᵢ * dᵢ +
            27 * bᵢ^2 * Vdiff +
            27 * aᵢ * dᵢ^2 +
            72 * aᵢ * cᵢ * Vdiff,
        )
        Q = ((Δ₁ + sqrt(Δ₁^2 - 4 * Δ₀^3)) / 2)^(1 / 3)
        S = sqrt(-2 * pᵢ / 3 + (Q + Δ₀ / Q) / (3 * aᵢ)) / 2

        root = sqrt(-4 * S^2 - 2 * pᵢ - qᵢ / S)

        # Check the 4 possible roots for being valid;
        # real and in [0,1]
        s1 = -bᵢ / (4 * aᵢ) + S + root / 2
        if valid(s1)
            return U_s(A, real(s1), i)
        end

        s2 = -bᵢ / (4 * aᵢ) + S - root / 2
        if valid(s2)
            return U_s(A, real(s2), i)
        end

        s3 = -bᵢ / (4 * aᵢ) - S + root / 2
        if valid(s3)
            return U_s(A, real(s3), i)
        end

        s4 = -bᵢ / (4 * aᵢ) - S - root / 2
        if valid(s4)
            return U_s(A, real(s4), i)
        end

        error("No valid root found, got $([s1,s2,s3,s4]).")
    else
        Vdiff = (V - A.t[idx - 1])
        @assert Vdiff >= 0

        # TODO: Consider degenerate case
        # where Δu[I+1] = 0
        i = Int((idx - 1) // 2)
        Δuᵢ₊₁ = Δu[i + 1]
        Δtᵢ₊₁ = Δt[i + 1]
        u_frac = u[i] / Δuᵢ₊₁
        root = sqrt(u_frac^2 + λ * (u_frac + λ / 4) + 2 * Vdiff / (Δtᵢ₊₁ * Δuᵢ₊₁))
        t[i] + (-u_frac + sign(u_frac) * root) * Δtᵢ₊₁
    end
end
