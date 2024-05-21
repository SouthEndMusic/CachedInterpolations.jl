`SmoothInterpolation.jl` exports 2 interpolation types in the style of [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl):

- `SmoothedLinearInterpolation`, a type of linear interpolation with well-behaved smoothed corners;
- `SmoothedLinearInterpolationIntInv`, the inverse of the antiderivative of a `SmoothedLinearInterpolation` if it exists.
