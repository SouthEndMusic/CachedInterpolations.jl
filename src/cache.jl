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

struct SmoothedLinearInterpolationIntInvCache{uType}
    # The degree of the polynomial whose roots need to be found
    degree::Vector{UInt8}
    # Quartic polynomial coefficients
    a::uType
    b::uType
    c::uType
    d::uType
    # Coefficients of depressed quartic
    p::uType
    q::uType
end

function SmoothedLinearInterpolationIntInvCache(A)
    coeffs = hcat(
        [
            collect(get_quartic_coefficients(A, idx)) for
            idx in eachindex(A.t) if idx ∉ [1, length(A.t)]
        ]...,
    )
    a, b, c, d = collect.(eachrow(coeffs))
    degree = fill(UInt8(4), length(a))

    p = @. (8 * a * c - 3 * b^2) / (8 * a^2)
    q = @. (b^3 - 4 * a * b * c + 8 * a^2 * d) / (8 * a^3)

    return SmoothedLinearInterpolationIntInvCache(degree, a, b, c, d, p, q)
end
