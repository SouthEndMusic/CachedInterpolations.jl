module SmoothInterpolation

using DataInterpolations:
    DataInterpolations, LinearInterpolation, AbstractInterpolation, munge_data, _interpolate
using FindFirstFunctions: searchsortedfirstcorrelated
using PrettyTables

include("cache.jl")
include("smoothed_linear_interpolation.jl")
include("integration_inverse.jl")
include("integration.jl")
include("utils.jl")

export LinearInterpolationIntInv,
    SmoothedLinearInterpolation,
    SmoothedLinearInterpolationIntInv,
    LinearInterpolation,
    invert_integral

end # SmoothInterpolation
