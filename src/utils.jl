"""
    S(A, t, idx)

Compute the spline parameter `s` from from the time `t`.

    ## Arguments

    - `A`: The `SmoothedLinearInterpolation` object
    - `t`: The time point
    - `idx`: The index indicating which spline section
"""
function S(A::SmoothedLinearInterpolation, t, idx)
    (; Δt, ΔΔt, degenerate_ΔΔt, t_tilde, λ) = A.cache
    Δtᵢ = Δt[idx]
    ΔΔtᵢ = ΔΔt[idx]
    tdiff = t - t_tilde[2 * idx - 1]
    @assert tdiff >= 0

    if degenerate_ΔΔt[idx]
        # Degenerate case Δtᵢ₊₁ ≈ Δtᵢ
        s = 1 / λ * tdiff / Δtᵢ
    else
        s = (-Δtᵢ + sqrt(Δtᵢ^2 + 2 * ΔΔtᵢ * tdiff / λ)) / ΔΔtᵢ
    end
    ε = 1e-5
    @assert -ε <= s <= 1 + ε "s should be in [0,1], got $s."
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

"""
Compute the coefficients for the quartic polynomial in s for the integration
of a spline section

Vdiff = c4 * s^4 + c3 * s^3 + c2 * s^2 + c1 * s + c0
"""
function get_quartic_coefficients(A::SmoothedLinearInterpolation, idx::Number)
    (; Δu, Δt, ΔΔu, ΔΔt, u_tilde, λ) = A.cache

    i = 2 * idx
    Δtᵢ = Δt[idx]
    Δuᵢ = Δu[idx]
    ΔΔuᵢ = ΔΔu[idx]
    ΔΔtᵢ = ΔΔt[idx]
    f₁ = u_tilde[i - 1] / λ
    f₂ = λ^2 / 24

    c4 = f₂ * (3 * ΔΔtᵢ * ΔΔuᵢ)
    c3 = f₂ * (4 * Δtᵢ * ΔΔuᵢ + 8 * ΔΔtᵢ * Δuᵢ)
    c2 = f₂ * (12 * ΔΔtᵢ * f₁ + 12 * Δtᵢ * Δuᵢ)
    c1 = f₂ * (24 * Δtᵢ * f₁)

    return c4, c3, c2, c1
end

"""
Determine whether a value s is valid, i.e.
- Its imaginary part is close to 0;
- Its real part is in the interval [0.1].
"""
valid(s) = isapprox(imag(s), 0; atol = 1e-4) && (0 <= real(s) <= 1)

"""
Compute t of a spline section from s.
"""
function T_s(A::SmoothedLinearInterpolationIntInv, s, idx)
    (; Δt, ΔΔt, t_tilde, λ) = A.cache
    Δtᵢ = Δt[idx]
    ΔΔtᵢ = ΔΔt[idx]
    return λ / 2 * ΔΔtᵢ * s^2 + λ * Δtᵢ * s + t_tilde[2 * idx - 1]
end

struct RootIterator{T1, T2, T3, T4, D}
    # T1: Complex constants
    # T2: Complex constants depending on c0,
    # so could carry differentiation information
    degree::D
    ab_part::T4
    c2::T2
    c3::T2
    Δ₀::T1
    Δ₁::T1
    Q::T1
    S::T1
    p::T3
    q::T3
end

"""
    iterate_roots(degree, c4, c3, c2, c1, c0, p, q)

Generate an iterator object which iterates over the roots
of the polynomial with the given coefficients of the given degree.
Coefficients for terms higher than the degree are not used, 
and p, q are only used when degree = 4.
"""
function iterate_roots(
    degree,
    c4::T2,
    c3::T2,
    c2::T2,
    c1::T2,
    c0::T1,
    p::T3,
    q::T3,
)::RootIterator where {T1, T2, T3}
    Δ₀ = zero(Complex(c0))
    Δ₁ = zero(Complex(c0))
    Q = zero(Complex(c0))
    S = zero(Complex(c0))
    if degree == 1
        ab_part = -c0 / c1
    elseif degree == 2
        Δ₀ = Complex(c1^2 - 4 * c2 * c0)
        ab_part = -c1 / (2 * c2)
    else
        Δ₀ = Complex(c2^2 - 3 * c3 * c1 + 12 * c4 * c0)
        Δ₁ = Complex(
            2 * c2^3 - 9 * c3 * c2 * c1 + 27 * c3^2 * c0 + 27 * c4 * c1^2 -
            72 * c4 * c2 * c0,
        )
        Q = ((Δ₁ + sqrt(Δ₁^2 - 4 * Δ₀^3)) / 2)^(1 / 3)
        if degree == 3
            ab_part = -c2 / (3 * c3)
        else
            ab_part = -c3 / (4 * c4)
            S = sqrt(-2 * p / 3 + (Q + Δ₀ / Q) / (3 * c4)) / 2
        end
    end
    return RootIterator(
        degree,
        complex(ab_part),
        complex(c2),
        complex(c3),
        Δ₀,
        Δ₁,
        Q,
        S,
        Complex(p),
        Complex(q),
    )
end

"""
Compute a root of a quartic polynomial
"""
function quartic_root(root_iterator::RootIterator{T1}, state)::T1 where {T1}
    (; ab_part, S, p, q) = root_iterator
    sign_1 = state < 3 ? -1 : 1
    sign_2 = (-1)^state
    root = sqrt(-4S^2 - 2p - sign_1 * q / S)
    out = ab_part + sign_1 * S + sign_2 * 0.5 * root
    return out
end

"""
Compute a root of a cubic polynomial
"""
function cube_root(root_iterator::RootIterator{T1}, state)::T1 where {T1}
    (; c3, Q, Δ₀, ab_part) = root_iterator
    ξ = exp(2π * im / 3)
    C = Q * ξ^(state - 1)
    return ab_part - (C + Δ₀ / C) / (3 * c3)
end

"""
Compute a root of a quadratic polynomial
"""
function square_root(root_iterator::RootIterator{T1}, state)::T1 where {T1}
    (; c2, ab_part, Δ₀) = root_iterator
    return ab_part + (-1)^state * sqrt(Δ₀) / (2 * c2)
end

"""
Compute a root of a linear polynomial
"""
function linear_root(root_iterator::RootIterator{T1}, state)::T1 where {T1}
    return root_iterator.ab_part
end

function root(root_iterator::RootIterator{T1}, state)::T1 where {T1}
    (; degree) = root_iterator
    if degree == 4
        quartic_root(root_iterator, state)
    elseif degree == 3
        cube_root(root_iterator, state)
    elseif degree == 2
        square_root(root_iterator, state)
    else
        linear_root(root_iterator, state)
    end
end

Base.length(root_iterator::RootIterator) = root_iterator.degree
Base.iterate(root_iterator::RootIterator) = (root(root_iterator, 1), 2)
Base.iterate(root_iterator::RootIterator, state) =
    state > root_iterator.degree ? nothing : (root(root_iterator, state), state + 1)

function p_coeff(c4::T2, c3::T2, c2::T2)::T2 where {T2}
    return (8 * c4 * c2 - 3 * c3^2) / (8 * c4^2)
end

function q_coeff(c4::T2, c3::T2, c2::T2, c1::T2)::T2 where {T2}
    return (c3^3 - 4 * c4 * c3 * c2 + 8 * c4^2 * c1) / (8 * c4^3)
end

"""
Compute u_tilde, the value of u at the boundary points between linear and spline sections
of a SmootedLinearInterpolation curve.
"""
function get_spline_ends(u, Δu, λ)
    u_tilde = zeros(2 * length(u))
    u_tilde[1] = u[1]
    u_tilde[2:2:(end - 1)] = u[1:(end - 1)] .+ (λ / 2) .* Δu[2:(end - 1)]
    u_tilde[3:2:end] = u[2:end] .- (λ / 2) .* Δu[2:(end - 1)]
    u_tilde[end] = u[end]
    return u_tilde
end