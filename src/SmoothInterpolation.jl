module SmoothInterpolation

using DataInterpolations:
    DataInterpolations, AbstractInterpolation, munge_data, _interpolate
using FindFirstFunctions: searchsortedfirstcorrelated

include("cache.jl")
include("smoothed_linear_interpolation.jl")
include("integration_inverse.jl")

export SmoothedLinearInterpolation, SmoothedLinearInterpolationIntInv

end # SmoothInterpolation
