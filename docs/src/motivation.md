# Motivation

I didn't like the options available in `DataInterpolations.jl` for my application, so I came up with my own concept that mainly keeps the linear interpolation intact, but rounds of the corners between the linear sections. The main recipe is this, per corner:
- Add 2 points on either side close to the corner, on their respectilve linear sections;
- Remove the corner point;
- Connect the 2 new points with a spline curve.

The advantage of the spline curve over a polynomial one is that the connection can be $C^1$ smooth (i.e. the smoothed curve and its derivative are continuous) without the possibility of introducing large oscillations with the use of degree 3 polynomials.