struct SmoothedLinearInterpolationIntInv{uType,tType,λType<:Real,T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    cache::SmoothedLinearInterpolationCache{uType,tType,λType}
    extrapolate::Bool
    function SmoothedLinearInterpolationIntInv(u, t, cache, λ, extrapolate)
        return new{typeof(u),typeof(t),typeof(λ),eltype(u)}(u, t, cache, extrapolate)
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
                0.5 * (u_tilde[end] + u_tilde[end-1]) * (t_tilde[end] - t_tilde[end-1])
        elseif j % 2 == 0
            i = Int(j // 2)
            Δuᵢ = Δu[i]
            Δuᵢ₊₁ = Δu[i+1]
            Δtᵢ = Δt[i]
            Δtᵢ₊₁ = Δt[i+1]
            U[j] =
                λ^2 / 24 * (Δtᵢ * (-3 * Δuᵢ + Δuᵢ₊₁) + Δtᵢ₊₁ * (-Δuᵢ + 3 * Δuᵢ₊₁)) +
                λ / 2 * (Δtᵢ + Δtᵢ₊₁) * u[i]
        else
            U[j] = 0.5 * (u_tilde[j] + u_tilde[j-1]) * (t_tilde[j] - t_tilde[j-1])
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
    t::Number,
    iguess,
)

    n_points = length(A.t)

    # idx of smallest idx such that A.t[idx] >= t
    idx = searchsortedfirstcorrelated(A.t, t, iguess)

    if idx == 1
        return 0
    elseif idx == 2

    elseif idx == n_points + 1

    elseif idx % 2 == 0

    else

    end
end
