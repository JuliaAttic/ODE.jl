using ODE
using ODE.ODETests
using Base.Test

@testset "ODE tests" begin
    include("iterators.jl")
    include("top-interface.jl")
    include("common.jl")
end
