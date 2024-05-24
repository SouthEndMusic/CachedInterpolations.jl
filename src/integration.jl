
function integrate_spline_section(A::SmoothedLinearInterpolation, idx::Number, t::Number)
    (; Δu, Δt, ΔΔu, ΔΔt, u_tilde, λ) = A.cache
    s = S(A, t, idx)

    i = 2 * idx
    Δtᵢ = Δt[idx]
    Δuᵢ = Δu[idx]
    ΔΔuᵢ = ΔΔu[idx]
    ΔΔtᵢ = ΔΔt[idx]
    f = u_tilde[i - 1] / λ

    a = 3 * ΔΔtᵢ * ΔΔuᵢ
    b = 4 * Δtᵢ * ΔΔuᵢ + 8 * ΔΔtᵢ * Δuᵢ
    c = 12 * ΔΔtᵢ * f + 12 * Δtᵢ * Δuᵢ
    d = 24 * Δtᵢ * f

    return λ^2 * (a * s^4 + b * s^3 + c * s^2 + d * s) / 24
end

DataInterpolations.samples(A::SmoothedLinearInterpolation) = (-1, 0)
function DataInterpolations._integral(
    A::SmoothedLinearInterpolation,
    idx::Number,
    t::Number,
)
    (; u_tilde, t_tilde) = A.cache

    if t == A.t[idx]
        return zero(eltype(A.u))
    end

    # idx of smallest idx such that A.t[idx] >= t
    idx = searchsortedfirstcorrelated(A.t, t, idx)

    i = 2 * idx
    u_tildeᵢ₋₃ = u_tilde[i - 3]
    u_tildeᵢ₋₂ = u_tilde[i - 2]

    t_tildeᵢ₋₃ = t_tilde[i - 3]
    t_tildeᵢ₋₂ = t_tilde[i - 2]

    # Integration of lower spline section
    if idx == 2
        # Special case of the first (half) spline section
        # which is linear
        if t <= t_tildeᵢ₋₂
            u_t = A(t)
            out = 0.5 * (t - t_tildeᵢ₋₃) * (u_t + u_tildeᵢ₋₃)
            return out
        else
            out = 0.5 * (t_tildeᵢ₋₂ - t_tildeᵢ₋₃) * (u_tildeᵢ₋₂ + u_tildeᵢ₋₃)
        end
    elseif idx == length(A.t) + 1
        # Special case of upper extrapolation
        u_t = A(t)
        out = 0.5 * (t - A.t[end]) * (u_t + A.u[end])
        return out
    else
        if t <= t_tildeᵢ₋₂
            out = integrate_spline_section(A, idx - 1, t)
            out -= integrate_spline_section(A, idx - 1, A.t[idx - 1])
            return out
        else
            out = integrate_spline_section(A, idx - 1, t_tildeᵢ₋₂)
            out -= integrate_spline_section(A, idx - 1, A.t[idx - 1])
        end
    end

    u_tildeᵢ₋₁ = u_tilde[i - 1]
    t_tildeᵢ₋₁ = t_tilde[i - 1]

    # Integration of linear section
    if t <= t_tildeᵢ₋₁
        u_t = A(t)
        out += 0.5 * (t - t_tildeᵢ₋₂) * (u_t + u_tildeᵢ₋₂)
        return out
    else
        out += 0.5 * (t_tildeᵢ₋₁ - t_tildeᵢ₋₂) * (u_tildeᵢ₋₁ + u_tildeᵢ₋₂)
    end

    # Integration of upper spline section
    out += integrate_spline_section(A, idx, t)
    return out
end
