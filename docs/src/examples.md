# Examples

## Smoothed linear interpolation

### An example of smooth linear interpolation

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

### The effect of the parameter λ

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