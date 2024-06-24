 [![codecov](https://codecov.io/gh/SouthEndMusic/CachedInterpolations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SouthEndmusic/CachedInterpolations.jl)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://southendmusic.github.io/CachedInterpolations.jl/dev/)

<img src="docs/src/assets/logo.svg" width="200">

# CachedInterpolations.jl

`CachedInterpolations.jl` exports interpolation types in the style of [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl):

- `CLinearInterpolation`, regular old linear interpolation;
- `CSmoothedLinearInterpolation`, a type of linear interpolation with well-behaved smoothed corners;
- `CSmoothedLinearInterpolationIntInv`, the inverse of the antiderivative of a `CSmoothedLinearInterpolation` if it exists;
- `CLinearInterpolationInvInv`, the inverse of the antiderivative of a `(C)LinearInterpolation` if it exists.

These interpolations make as much use as possible of cached data computed at initialization, for a speed against memory usage tradeoff.

## Installation

Simply run the code below from the package manager.

```
pkg> add CachedInterpolations
```

## Supported features

Not all features for interpolation objects from `DataInterpolations.jl` are currently supported.

|                                      | Evaluation | Derivative    | Integration                                |
| -----------------------------------  | ---------- | ------------- | ------------------------------------------ |
| `CLinearInterpolation`               | Supported  | Not supported | Supported                                  |
| `CSmoothedLinearInterpolation`       | Supported  | supported     | Supported                                  |
| `CSmoothedLinearInterpolationIntInv` | Supported  | supported     | Not supported                              |
| `CLinearInterpolationIntInv`         | Supported  | supported     | Not supported                              |

If you wish to use one of the currently unsupported features, please [write an issue](https://github.com/SouthEndMusic/CachedInterpolations.jl/issues). Note that differentiation can also be achieved with many of the [automatic differentiation packages](https://juliadiff.org/#the_big_list) in the Julia ecosystem.

## Logo

The logo is inspired by the [julia logo graphics](https://github.com/JuliaLang/julia-logo-graphics).
