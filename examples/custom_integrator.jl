"""

Here we demonstrate how to implement an integrator so that it is
compatible with `ODE.jl`.  The implementation can be contained in an
external package, but if the API conforms to our backend one could use
this integrator to solve an explicit ODE with `solve(...)`.

"""

module MyIntegrator

using ODE
# We have to define these methods on our integrator
import ODE: init, onestep!
import ODE: AbstractIntegrator, AbstractState, ExplicitODE, Step

# The options are stored in the integrator type.  Integrator is
# immutable: we can't change the options once we start the
# integration.
immutable EulerIntegrator{T} <: AbstractIntegrator{T}
    tstop::T
    initstep::T
end

# This is necessary for the `solve` to construct the iterator.  Every
# iterator has to have a constructor of the form
# `Iterator(::IVP;opts...)`, here we implement an interator that only
# works with explicit differential equations, hence we restrict the
# constructor to `ExplicitODE`.  If someone tries to use our iterator
# to solve an equation of unsupported type `solve` will throw an error
# (there is a fallback constructor for `AbstractIntegrator` that
# throws an error).
function EulerIntegrator{T}(ode::ExplicitODE{T};
                            tstop = T(Inf),
                            initstep = T(1//10),
                            opts...)
    EulerIntegrator{T}(tstop,initstep)
end

# The state of the integrator, it stores the current values for time
# `t` and the solution `y` itself, along with its derivative,
# `dy`.  For convenience we used the already available type to store
# this kind of information, `Step`, that has three fields: `t`, `y`
# and `dy`.  See the end of the module for an alternative definition
# of a state.
type EulerState{T,Y} <: AbstractState{T,Y}
    step::ODE.Step{T,Y}
end
output(state::EulerState) = output(state.step)

# Generates the state given an ODE and an integrator, this method has
# to be specialized for `EulerIntegrator`, we don't have to specialize
# it for `ExplicitODE` as the type of an ODE is already filtered by
# the specialized the construtor, but we do it here for clarity.
function init(ode::ExplicitODE, integ::EulerIntegrator)
    EulerState(ODE.Step(ode.t0,copy(ode.y0),copy(ode.dy0)))
end

function onestep!(ode::ExplicitODE, integ::EulerIntegrator, state::EulerState)
    # as mentioned before, this function unpacks the variables
    # `t,y,dy` from the current state.  It is not necessary, you could
    # access the fields directly but we use it here for convenience.
    t, y, dy = ODE.output(state)

    tdir = sign(integ.tstop-ode.t0)

    # the only stop condition our solver has.  Note the use of `abs`,
    # which enables integration backward in time.
    if tdir*t >= tdir*integ.tstop
        # this flag finalizes the iterator
        return ODE.finish
    else
        # trim the stepsize to match the `tstop`, prevents
        # overshooting
        dt  = tdir*min(integ.initstep,abs(integ.tstop-t))

        # update the time,
        state.step.t += dt
        # the function (this is the `dy` from the previous step)
        state.step.y += dt*dy
        # and its derivative
        ode.F!(t,y,dy)
        # return a flag to continue the integration
        return ODE.cont
    end
end

# OPTIONAL: Define properties of this integrator: order, name and
# whether it is adaptive or not.  At this point the information
# supplied here is not used but it might be a good idea to implement
# these methods for future use.
order{T}(::Type{EulerIntegrator{T}}) = 1
name{T}(::Type{EulerIntegrator{T}}) = "My own Euler integrator"
isadaptive{T}(::Type{EulerIntegrator{T}}) = false

# OPTIONAL:
# Another possiblity to implement state would be to declare
type EulerState2{T,Y} <: ODE.AbstractState{T,Y}
    t::T
    y::Y
    dy::Y
end
# but then we have to define the method `output`
output(state::EulerState2) = (state.t, state.y, state.dy)

end

# Usage example
using ODE
using ODE.ODETests
using MyIntegrator

integ = MyIntegrator.EulerIntegrator

# declare the ODE as usual
ode = ODE.ExplicitODE(0.0,[1.0],(t,y,dy)->copy!(dy,y))
# solve the `ode` with our integrator, note that we can pass options to `solve`
sol = ODE.solve(ode,solver=integ, tstop=1.0, initstep=0.001)
# print the last step of the solution

# the solution can be accessed with
(sol.t, sol.y, sol.dy)

# test the integrator
# TODO: for now I can't figure out why it fails.
ODETests.test_integrator(integ)
