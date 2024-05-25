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
        s = 1 / λ * tdiff / Δtᵢ
    else
        s = (-Δtᵢ + sqrt(Δtᵢ^2 + 2 * ΔΔtᵢ * tdiff / λ)) / ΔΔtᵢ
    end
    @assert 0 <= s <= 1 "s should be in [0,1], got $s."
    return s
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
