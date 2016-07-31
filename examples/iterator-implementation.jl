using ODE

# our test ODE:
t0 = 0.0
y0 = [1.0]
ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])

# let's begin

type EulerIntegrator <: ODE.AbstractIntegrator
    # nothing in here yet
end
# EulerIntegrator works with ExplicitODE only:
EulerIntegrator(ode::ODE.ExplicitODE; opts...) = EulerIntegrator()

type EulerState
    t
    y
    dy
end

ODE.output(state::EulerState) = state.t, state.y, state.dy

function ODE.init(ode::ODE.ExplicitODE, integ::EulerIntegrator)
    t0, y0 = ode.t0, ode.y0
    dy0 = similar(ode.dy0)
    ode.F!(t0,y0,dy0)           # fill in the values of the derivative
    EulerState(t0,y0,dy0)
end

function ODE.onestep!(ode::ODE.ExplicitODE, integ::EulerIntegrator, state::EulerState)
    t, y, dy = ODE.output(state)
    dt = 0.1
    state.y += dt*dy
    state.t += dt
    ode.F!(t,y,dy)           # update the derivative
    return ODE.cont
end

# use it:
prob = ODE.Problem(ode, EulerIntegrator)
for (t,y) in prob
    if t > 1
        println("$(prob.solver) $t $y")
        break
    end
end



"""

There are several problems with the above implementation.  First of
all, it has a constant prescribed step size.  This could easily be
fixed by changing the type definition and the `solver` to

"""

type EulerIntegrator2 <: ODE.AbstractIntegrator
    initstep
end
EulerIntegrator2(ode::ODE.ExplicitODE; initstep=0.1, opts...) = EulerIntegrator2(initstep)

"""
we should also change the line `dt = 0.1` in the `onestep!` function
to `dt = stepper.initstep`.  Now we can run our integrator with a
custom step size!
"""

function ODE.init(ode::ODE.ExplicitODE, integ::EulerIntegrator2)
    t0, y0 = ode.t0, ode.y0
    dy0 = similar(ode.dy0)
    ode.F!(t0,y0,dy0)           # fill in the values of the derivative
    EulerState(t0,y0,dy0) # we re-use the state
end

function ODE.onestep!(ode::ODE.ExplicitODE, integ::EulerIntegrator2, state::EulerState)
    t, y, dy = ODE.output(state)
    dt = 0.1
    state.y += dt*dy
    state.t += dt
    ode.F!(t,y,dy)           # update the derivative
    return ODE.cont
end

# use it:
prob = ODE.Problem(ode, EulerIntegrator2)
for (t,y) in prob
    if t > 1
        println("$(prob.solver) $t $y")
        break
    end
end



"""

Another issue is type stability, to make `EulerIntegrator` perform
better we should explicitly annotate the fields in both
`EulerIntegrator` and `EulerState` like this

"""

type EulerIntegrator3{T} <: ODE.AbstractIntegrator{T}
    initstep::T
end
EulerIntegrator3{T}(ode::ODE.ExplicitODE{T}; initstep::T = T(0.1), opts...) = EulerIntegrator3(initstep)

type EulerState3{T,Y}
    t::T
    y::Y
    dy::Y
end

"""

But the `EulerState{T,Y}` is exactly the same as `Step` from
`base.jl`, so we can simplify it a bit more

"""

type EulerState33{T,Y}
    step::ODE.Step{T,Y}
end
ODE.output(state::EulerState33) = ODE.output(state.step)

function ODE.init(ode::ODE.ExplicitODE, integ::EulerIntegrator3)
    t0, y0 = ode.t0, ode.y0
    dy0 = similar(ode.dy0)
    ode.F!(t0,y0,dy0)           # fill in the values of the derivative
    EulerState33(ODE.Step(t0,y0,dy0)) # we re-use the state
end

function ODE.onestep!(ode::ODE.ExplicitODE, integ::EulerIntegrator3, state::EulerState33)
    t, y, dy = ODE.output(state)
    dt = 0.1
    state.step.y += dt*dy
    state.step.t += dt
    ode.F!(t,y,dy)           # update the derivative
    return ODE.cont
end

# use it:
prob = ODE.Problem(ode, EulerIntegrator3)
for (t,y) in prob
    if t > 1
        println("$(prob.solver) $t $y")
        break
    end
end




"""

Once we do that, in the `init` we should change
`EulerState(t0,y0,dy0)` to `EulerState(Step(t0,y0,dy0))` and redefine
`output` to `output(state::EulerState)=output(state.step)`
(`output(::Step)` is already defined in `base.jl`).

One could even replace `EulerState` with `Step` completely, but this
won't allow us to extend the state with some additional variables and
storage space in the future.

The last thing is that our stepper will continue the integration
forever: it doesn't have a stopping condition.  We could add one as an
option.

"""

type EulerIntegrator4{T} <: ODE.AbstractIntegrator{T}
    initstep::T
    tstop::T
end
EulerIntegrator4{T}(ode::ODE.ExplicitODE{T};
                    initstep::T = T(0.1),
                    tstop::T = T(Inf),
                    opts...) = EulerIntegrator4(initstep,tstop)

function ODE.init(ode::ODE.ExplicitODE, integ::EulerIntegrator4)
    t0, y0 = ode.t0, ode.y0
    dy0 = similar(ode.dy0)
    ode.F!(t0,y0,dy0)           # fill in the values of the derivative
    EulerState33(ODE.Step(t0,y0,dy0)) # we re-use the state
end

function ODE.onestep!(ode::ODE.ExplicitODE, integ::EulerIntegrator4, state::EulerState33)
    t, y, dy = ODE.output(state)

    if t > integ.tstop
        return ODE.finished
    end

    dt = integ.initstep
    state.step.y += dt*dy
    state.step.t += dt
    ode.F!(t,y,dy)           # update the derivative
    return ODE.cont
end

# use it:
prob = ODE.Problem(ode, EulerIntegrator4)
for (t,y) in prob
    if t > 1
        println("$(prob.solver) $t $y")
        break
    end
end

"""

As a final improvement, we can (although this is not necessary) use a
structure `FixedOptions` from `options.jl` to keep our options in one
structure.  A corresponding options type for adaptive solver is
`AdaptiveOptions`.  This way we can use the standarized defaults for
most options and keep our solver in line with the standard
keywords.  Naturally, we have to update `onestep!` to use the subtype.

"""

type EulerIntegrator5{T} <: ODE.AbstractIntegrator{T}
    opts::ODE.FixedOptions{T}
    tstop::T
end
function EulerIntegrator5{T}(ode::ODE.ExplicitODE{T};
                             tstop::T = T(Inf),
                             opts...)
    opts_fixed = ODE.FixedOptions{T}(;opts...)
    EulerIntegrator5(opts_fixed, tstop)
end

function ODE.init(ode::ODE.ExplicitODE, integ::EulerIntegrator5)
    t0, y0 = ode.t0, ode.y0
    dy0 = similar(ode.dy0)
    ode.F!(t0,y0,dy0)           # fill in the values of the derivative
    EulerState33(ODE.Step(t0,y0,dy0)) # we re-use the state
end

function ODE.onestep!(ode::ODE.ExplicitODE, integ::EulerIntegrator5, state::EulerState33)
    t, y, dy = ODE.output(state)

    if t > integ.tstop
        return ODE.finished
    end

    dt = integ.opts.initstep
    state.step.y += dt*dy
    state.step.t += dt
    ode.F!(t,y,dy)           # update the derivative
    return ODE.cont
end

# use it:
prob = ODE.Problem(ode, EulerIntegrator5)
for (t,y) in prob
    if t > 1
        println("$(prob.solver) $t $y")
        break
    end
end
