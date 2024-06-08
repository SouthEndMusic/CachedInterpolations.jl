# Examples

## Smoothed linear interpolation

```@example 1
import Random # hide
Random.seed!(2) # hide
using SmoothInterpolation

u = rand(10)
t = cumsum(rand(10))

itp = SmoothedLinearInterpolation(u, t; extrapolate = true)
```

```@example 1
using Plots
scatter(itp.t, itp.u, label = "Input")
plot!(itp, label = "Smoothed Linear Interpolation")
```

## Inverting the integral

```@example 1
itp_int_inv = invert_integral(itp)
V = 1.0
t_V = itp_int_inv(V)
t_eval_V = range(t[1], t_V, length = 100)
plot!(t_eval_V, itp.(t_eval_V), fill = (:blue, 0, 0.5), label = "area = $V")
```

!!! tip
    The integral inverse of `SmoothedLinearInterpolation` is expensive to compute as it involves solving a quartic equation. If performance is important to your application, consider converting your `SmoothedLinearInterpolation` object into a `LinearInterpolation` object using `LinearInterpolation(A::SmoothedLinearInterpolation; n_samples = 10)`, which samples the spline sections. The inverse of this is much cheaper.


## The effect of the parameter λ

```@example 1
using ColorSchemes
t = [0, 1, 2, 2.5, 3, 3.5, 4]
u = Float64[-1, 1, -1, 0, 1, 0, -1]
pl = plot()
scatter!(t, u, label = "Input", legend = :top)
N = 101
Λ = range(0, 1, length = N)
colors = cgrad(:jet, range(0, 1, length = N))

for (i, (λ, color)) in enumerate(zip(Λ, colors))
    itp = SmoothedLinearInterpolation(u, t; λ)
    label = i % 10 == 1 ? "λ = $λ" : nothing 
    plot!(itp; label, color)
end

pl
```

## Derivatives

Derivatives can be calculated using `DataInterpolations.derivative(itp, t)`. There is a quite simple relationship between the derivative of the inverse of the integral of a function and the function itself:

```math
(F^{-1})'(V) = \frac{1}{f(F^{-1}(V))}.
```


See also the code example below.

```@example 1
using DataInterpolations
using ForwardDiff
Random.seed!(15) # hide

t = cumsum(rand(10))
u = rand(10)

itp = SmoothedLinearInterpolation(u, t; extrapolate = true)
itp_int_inv = invert_integral(itp)
u_int_eval = itp_int_inv.t[1]:0.01:(itp_int_inv.t[end] + 1)

# Compute the hardcoded SmoothedLinearInterpolationIntInv derivative
t_deriv_eval = DataInterpolations.derivative.(Ref(itp_int_inv), u_int_eval)

# Compute the SmoothedLinearInterpolationIntInv derivative using ForwardDiff
t_deriv_forward_diff = ForwardDiff.derivative.(Ref(itp_int_inv), u_int_eval)

# Compare results
@show t_deriv_eval ≈ 1 ./ itp.(itp_int_inv.(u_int_eval))
@show t_deriv_eval ≈ t_deriv_forward_diff
```