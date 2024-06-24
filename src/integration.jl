"""
    integrate_spline_section(A::CSmoothedLinearInterpolation, idx::Number, t::Number)

Integrate the idx-th spline section from its lower time bound up to t

## Arguments

 - `A`: CSmoothedLinearInterpolation object
 - `idx`: Index of the spline section
 - `t`: upper integration bound
"""
function integrate_spline_section(A::CSmoothedLinearInterpolation, idx::Number, t::Number)
    s = S(A, t, idx)
    c4, c3, c2, c1 = get_quartic_coefficients(A, idx)

    return c4 * s^4 + c3 * s^3 + c2 * s^2 + c1 * s
end

function DataInterpolations._integral(
    A::CSmoothedLinearInterpolation,
    idx::Number,
    t::Number,
)
    (; u_tilde, t_tilde, idx_prev) = A.cache

    if t == A.t[idx]
        return zero(eltype(A.u))
    end

    # idx of smallest idx such that A.t[idx] >= t
    idx = searchsortedfirstcorrelated(A.t, t, idx_prev[])
    idx_prev[] = idx

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
        out = 0.5 * (t - A.t[end - 1]) * (u_t + A.u[end - 1])
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

function DataInterpolations._integral(A::CLinearInterpolation, idx::Number, t::Number)
    (; slope) = A.cache
    idx += 1
    tdiff = (t - A.t[idx])
    return tdiff * (A.u[idx] + 0.5 * slope[idx] * tdiff)
end