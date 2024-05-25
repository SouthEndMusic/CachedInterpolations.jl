using SafeTestsets

@safetestset "Interpolation" include("interpolation.jl")
@safetestset "Extrapolation" include("extrapolation.jl")
@safetestset "Continuity" include("continuity.jl")
@safetestset "Integration" include("integration.jl")