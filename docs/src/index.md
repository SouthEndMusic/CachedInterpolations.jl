# SmoothInterpolation.jl

`SmoothInterpolation.jl` exports 2 interpolation types in the style of [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl):

- `CSmoothedLinearInterpolation`, a type of linear interpolation with well-behaved smoothed corners;
- `CSmoothedLinearInterpolationIntInv`, the inverse of the antiderivative of a `CSmoothedLinearInterpolation` if it exists;
- `CLinearInterpolationInvInv`, the inverse of the antiderivative of a `(C)LinearInterpolation` if it exists.

## Installation

Simply run the code below from the package manager.

```
pkg> add SmoothInterpolation
```

## Supported features

Not all features for interpolation objects from `DataInterpolations.jl` are currently supported.

|                                     | Evaluation | Derivative    | Integration                                |
| ----------------------------------- | ---------- | ------------- | ------------------------------------------ |
| `CSmoothedLinearInterpolation`       | Supported  | supported     | Supported                                  |
| `CSmoothedLinearInterpolationIntInv` | Supported  | supported     | Not supported                              |
| `CLinearInterpolationIntInv`         | Supported  | supported     | Not supported                              |

If you wish to use one of the currently unsupported features, please [write an issue](https://github.com/SouthEndMusic/SmoothInterpolation.jl/issues). Note that differentiation can also be achieved with many of the [automatic differentiation packages](https://juliadiff.org/#the_big_list) in the Julia ecosystem.

## Logo

The logo is inspired by the [julia logo graphics](https://github.com/JuliaLang/julia-logo-graphics).