using SafeTestsets

@safetestset "Interpolation" include("interpolation.jl")
@safetestset "Extrapolation" include("extrapolation.jl")