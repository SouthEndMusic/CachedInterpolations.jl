using SafeTestsets

@safetestset "Utils" include("utils_test.jl")
@safetestset "Interpolation" include("interpolation_test.jl")
@safetestset "Extrapolation" include("extrapolation_test.jl")
@safetestset "Continuity" include("continuity_test.jl")
@safetestset "Integration" include("integration_test.jl")
@safetestset "Derivatives" include("derivatives_test.jl")