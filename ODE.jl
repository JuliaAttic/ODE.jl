# This is a prototype for new ODE.jl internals.

using Compat
import Compat.view
using Parameters
using ForwardDiff
# ob: packages for multistep methods
using DataStructures
using Polynomials

include("helpers.jl")

# core.jl contains the core iteration interface
include("core.jl")
# the dense output
include("dense.jl")
include("events.jl")

# ui.jl contains all the utilities to make the core.jl usable.
include("ui.jl")

#
include("tableaus.jl")

###
# Solvers
# call into the API provided in core.jl
include("runge_kutta.jl")
include("runge_kutta_v2.jl")

include("rosenbrock.jl")
include("fixed_adam.jl")
