function DataInterpolations._derivative(
    A::SmoothedLinearInterpolation{<:AbstractVector},
    t::Number,
    iguess,
)
    (; u, t_tilde, idx_prev) = A.cache

    # idx of smallest idx such that A.t[idx] >= t
    idx = searchsortedfirstcorrelated(A.t, t, idx_prev[])
    A.cache.idx_prev[] = idx

    if idx == 1 || idx == length(u) + 1
        # Linear extrapolation
        A.cache.linear_slope[idx]
    else
        # Interpolation
        if t < t_tilde[2 * idx - 2]
            U_deriv(A, t, idx - 1)
        elseif t > t_tilde[2 * idx - 1]
            U_deriv(A, t, idx)
        else
            # Linear interpolation
            A.cache.linear_slope[idx]
        end
    end
end

function DataInterpolations._derivative(A::AbstractInterpolationIntInv, V::Number, iguess)
    t = A(V, iguess)
    return 1 / forward_itp(A)(t)
end
