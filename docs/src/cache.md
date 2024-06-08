# To Cache or not to Cache

At the initialization of the interpolation objects exposed by this package, a lot of data is precomputed and cached, see for instance the example below.

```@example 1
import Random # hide
Random.seed!(2) # hide
using SmoothInterpolation

u = rand(10)
t = cumsum(rand(10))

itp = SmoothedLinearInterpolation(u, t)
itp.cache
```

This means that evaluation of the interpolation is faster, at the cost of more memory allocation at initialization. This is in contrast to the interpolation in `DataInterpolations.jl`, where very little to no memory is allocated at the initialization of interpolation objects. What is better depends on the application. 

If you want to use the interpolation objects exposed by this package without pre-allocation, please let me know in [this issue](https://github.com/SouthEndMusic/SmoothInterpolation.jl/issues/45).