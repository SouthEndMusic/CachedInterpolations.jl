function integrate_spline_section(A::SmoothedLinearInterpolation, idx::Number, t::Number) end

function DataInterpolations._integral(
    A::SmoothedLinearInterpolation{uType},
    idx::Number,
    t::Number,
)
    (; u_tilde, t_tilde) = A.cache
    out = zero(uType)
    i = 2 * idx

    t_tildeᵢ₋₂ = t_tilde[i - 2]

    if t < t_tildeᵢ₋₂
        out += integrate_spline_section(A, ..., t_tildeᵢ₋₂)
        out -= integrate_spline_section(A, ..., A.t[idx])
        return out
    else
        t_tildeᵢ₋₁ = t_tilde[i - 1]
        out += integrate_spline_section(A, ..., t_tildeᵢ₋₂)
        out -= integrate_spline_section(A, ..., t_tildeᵢ₋₁)
    end

    t_tildeᵢ = t_tilde[i]

    if t < t_tildeᵢ
        u_t = A(t)
        u_tildeᵢ₋₁ = u_tilde[i - 1]
        out += 0.5 * (t - t_tildeᵢ₋₁) * (u_t + u_tildeᵢ₋₁)
        return out
    else
        u_tildeᵢ₋₁ = u_tilde[i]
        out += 0.5 * (t_tildeᵢ - t_tildeᵢ₋₁) * (u_tildeᵢ + u_tildeᵢ₋₁)
    end

    out += integrate_spline_section(A, ..., t)
    return out
end