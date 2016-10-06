using ODE
using ODE.ODETests
using Base.Test
using ForwardDiff

@testset "ODE tests" begin
    include("iterators.jl")
    include("top-interface.jl")
    include("dense.jl")
end
