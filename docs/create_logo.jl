using SmoothInterpolation
using Luxor

height = 500
wide = false
Drawing(wide ? 2 * height : height, height, normpath(@__DIR__, "src/assets/logo.svg"))
origin()

# Distance of dots from the origin
R = 0.4 * height

# Radius of dots
r = 0.05 * height

# Dot coordinates
phi = [(i + 1) * 2π / 3 - π / 2 for i in 1:3]
x_dots = R * cos.(phi)
y_dots = R * sin.(phi)

# Interpolation curves
N = 100
setline(10)
for λ in 0.5:0.25:1.0
    itp = CSmoothedLinearInterpolation(y_dots, x_dots; λ, extrapolate = true)
    x_curve = range(x_dots[1], x_dots[end]; length = N)
    y_curve = itp.(x_curve)
    points = Point.(x_curve, y_curve)
    setcolor(λ .* Luxor.julia_blue)
    for i in 1:(N - 1)
        line(points[i + 1], points[i], :stroke)
    end
end

# The julia colored dots
colors = [Luxor.julia_purple, Luxor.julia_red, Luxor.julia_green]
for (x, y, c) in zip(x_dots, y_dots, colors)
    setcolor(c)
    circle(Point(x, y), r; action = :fill)
end

finish()
preview()
