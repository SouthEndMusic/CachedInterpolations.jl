struct SmoothedLinearInterpolationCache{uType, tType, λType <: Number}
    u::uType
    t::tType
    Δu::uType
    Δt::tType
    ΔΔu::uType
    ΔΔt::tType
    u_tilde::uType
    t_tilde::tType
    linear_slope::uType
    λ::λType
end

function get_spline_ends(u, Δu, λ)
    u_tilde = zeros(2 * length(u))
    u_tilde[1] = u[1]
    u_tilde[2:2:(end - 1)] = u[1:(end - 1)] .+ (λ / 2) .* Δu[2:(end - 1)]
    u_tilde[3:2:end] = u[2:end] .- (λ / 2) .* Δu[2:(end - 1)]
    u_tilde[end] = u[end]
    return u_tilde
end

function SmoothedLinearInterpolationCache(u, t, λ)::SmoothedLinearInterpolationCache
    Δu = diff(u)
    Δt = diff(t)
    @assert !any(iszero.(Δt))
    pushfirst!(Δt, Δt[1])
    push!(Δt, Δt[end])
    pushfirst!(Δu, Δu[1])
    push!(Δu, Δu[end])
    ΔΔu = diff(Δu)
    ΔΔt = diff(Δt)
    u_tilde = get_spline_ends(u, Δu, λ)
    t_tilde = get_spline_ends(t, Δt, λ)
    linear_slope = Δu ./ Δt
    return SmoothedLinearInterpolationCache(
        u,
        t,
        Δu,
        Δt,
        ΔΔu,
        ΔΔt,
        u_tilde,
        t_tilde,
        linear_slope,
        λ,
    )
end
