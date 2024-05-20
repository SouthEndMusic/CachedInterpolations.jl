struct SmoothedLinearInterpolationIntInv{uType, tType, λType <: Real, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    cache::SmoothedLinearInterpolationCache{uType, tType, λType}
    extrapolate::Bool
    function SmoothedLinearInterpolationIntInv(u, t, cache, λ, extrapolate)
        return new{typeof(u), typeof(t), typeof(λ), eltype(u)}(u, t, cache, extrapolate)
    end
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

function SmoothedLinearInterpolationIntInv(
    itp::SmoothedLinearInterpolation,
)::SmoothedLinearInterpolationIntInv
    (; cache, extrapolate) = itp
    t = integrate_sections(cache)
    u = cache.t_tilde
    return SmoothedLinearInterpolationIntInv(u, t, cache, cache.λ, extrapolate)
end

function DataInterpolations._interpolate(
    A::SmoothedLinearInterpolationIntInv{<:AbstractVector},
    V::Number,
    iguess,
)
    n_points = length(A.t)
    (; u, t, Δu, Δt, ΔΔu, ΔΔt, u_tilde, λ) = A.cache

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
        Δtᵢ = Δt[i]
        Δuᵢ = Δu[i]
        ΔΔuᵢ = ΔΔu[i]
        ΔΔtᵢ = ΔΔt[i]

        f = Δuᵢ + u_tilde[idx - 1] / λ

        # Independent of V, can be cached
        a = 3 * ΔΔtᵢ * ΔΔuᵢ
        b = 4 * Δtᵢ * ΔΔuᵢ
        c = 12 * ΔΔtᵢ * f
        d = 24 * Δtᵢ * f

        e = complex(-24 * Vdiff / λ^2)

        # Independent of V, can be cached
        p = (8 * a * c - 3 * b^2) / (8 * a^2)
        q = (b^3 - 4 * a * b * c + 8 * a^2 * d) / (8 * a^3)

        Δ₀ = c^2 - 3 * b * d + 12 * a * e
        Δ₁ = 2 * c^3 - 9 * b * c * d + 27 * b^2 * e + 27 * a * d^2 - 72 * a * c * e
        Q = ((Δ₁ + sqrt(Δ₁^2 - 4 * Δ₀^3)) / 2)^(1 / 3)
        S = sqrt(-2 * p / 3 + (Q + Δ₀ / Q) / (3 * a)) / 2

        root = sqrt(-4 * S^2 - 2 * p - q / S)
        s = real(-b / (4 * a) + S + root / 2)
        if !(0 <= s <= 1)
            s = real(-b / (4 * a) + S - root / 2)
        end
        if !(0 <= s <= 1)
            s = real(-b / (4 * a) - S + root / 2)
        end
        if !(0 <= s <= 1)
            s = real(-b / (4 * a) - S - root / 2)
        end
        @assert (0 <= s <= 1)
        U_s(A, s, i)
        a * s^4 + b * s^3 + c * s^2 + d * s + e
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
