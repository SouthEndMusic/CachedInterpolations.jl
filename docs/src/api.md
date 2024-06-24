# API

Documentation for `CachedInterpolations.jl`'s public interface.

```@docs
CSmoothedLinearInterpolation
CLinearInterpolation(::CSmoothedLinearInterpolation)
invert_integral(::Union{LinearInterpolation, CLinearInterpolation})
invert_integral(::CSmoothedLinearInterpolation)
```