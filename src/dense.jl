# A higher level stepper, defined as a wrapper around another stepper.

"""

Dense output options:

- tspan    ::Vector{T}  output times
- points   ::Symbol which points are returned: `:specified` only the
  ones in tspan or `:all` which includes also the step-points of the solver.
- stopevent  Stop integration at a zero of this function
- roottol    TODO

"""

immutable DenseOptions{T<:Number} <: Options{T}
    tout::Vector{T}
    # points   ::Symbol
    # stopevent::S
    # roottol  ::T
end

@compat function (::Type{DenseOptions{T}}){T}(;
                                              tstop        = T(Inf),
                                              tout::Vector = T[tstop],
                                              # points::Symbol= :all,
                                              # stopevent::S  = (t,y)->false,
                                              # roottol       = eps(T)^T(1//3),
                                              kargs...)
    DenseOptions{T}(tout)
end


"""

A stepper specialized in dense output.  It wraps around another
`Problem` and stores the subsequent steps generated by `Problem` and
interpolates the results on request (currently this means at the
output times stored in `options.tout`).

"""
immutable DenseOutput{P<:Problem,OP<:DenseOptions} <: AbstractSolver
    prob::P
    options::OP
end

function solve{S<:DenseOutput}(ivp::IVP,
                               ::Type{S};
                               method = RKIntegratorAdaptive{:rk45},
                               options...)
    T = eltype(ivp)[1]
    sol_orig = Problem(ivp,method{T}(; options...))
    dense_options = DenseOptions{T}(; options...)
    dense_stepper = S(sol_orig,dense_options)
    return Problem(ivp,dense_stepper) # TODO: this is where it is needed.
end

"""

The state of the dense stepper

"""
type DenseState{St<:AbstractState,T,Y} <: AbstractState{T,Y}
    tout_i::Int
    step_prev::Step{T,Y}
    step_out::Step{T,Y}
    integrator_state::St
end

output(ds::DenseState) = output(ds.step_out)

function init(ivp::IVP,
              stepper::DenseOutput)
    ivp = stepper.prob.ivp
    integrator_state = init(stepper.prob.ivp, stepper.prob.stepper)
    dy0 = similar(ivp.y0)
    ivp.F!(ivp.t0,ivp.y0,dy0)
    step_prev = Step(ivp.t0,copy(ivp.y0),dy0)
    step_out = Step(ivp.t0,similar(ivp.y0),similar(ivp.y0))
    return DenseState(1,step_prev,step_out,integrator_state)
end


"""

TODO: rename `tout` to `tout` and drop the support for
`points=:all` outside of the `odeXX`?  Maybe even
`odeXX(;tout=[...])` would use dense output while `odeXX(;)`
wouldn't.

"""

function onestep!(ivp::IVP,
                  stepper::DenseOutput,
                  state::DenseState)
    i = state.tout_i
    if i > length(stepper.options.tout)
        return finish
    end

    # our next output time
    ti = stepper.options.tout[i]

    sol = stepper.prob # this looks weird
    sol_state = state.integrator_state

    # try to get a new set of steps enclosing `ti`, if all goes
    # right we end up with t∈[t1,t2] with
    # t1,_=output(state.step_prev)
    # t2,_=output(state.integrator_state)
    status = next_interval!(sol,sol_state,state.step_prev,ti)
    if status == abort
        # we failed to get enough steps
        warn("Iterator was exhausted before the dense output could produce the output.")
        return abort
    else
        # we got the steps, proceed with the interpolation, this fills
        # the state.step_out with y(ti) and y'(ti) according to an
        # interpolation algorithm specific for a method (defaults to
        # hermite O(3)).
        interpolate!(state.integrator_state,state.step_prev,ti,state.step_out)

        # increase the counter
        state.tout_i += 1
        return cont
    end
end

"""

Pulls the results from the (prob,state) pair using `onestep!` until
we reach a first step such that `t>=tout`.  It fills the `steps`
variable with (Step(t1,y(t1),dy(t1)),Step(t2,y(t2),dy(t2))), where
`t1` is is the step before `tout` and `t2` is `>=tout`.  In
other words `tout∈[t1,t2]`.

TODO: tdir

"""
function next_interval!(prob,state,step_prev,tout)

    while true
        # get the current time
        t1   = step_prev.t
        t2,_ = output(state)
        t1, t2 = sort([t1,t2])
        if t1 <= tout <= t2
            # we found the enclosing times
            return cont
        end

        # save the current state of solution
        t, y, dy = output(state)
        step_prev.t = t
        copy!(step_prev.y,y)
        copy!(step_prev.dy,dy)

        # try to perform a single step with the prob
        status = onestep!(prob.ivp, prob.stepper, state)

        if status != cont
            return status
        end
    end

    # this will never happen
    return abort
end

"""

Make dense output using Hermite interpolation of order O(3). Updates
yout in-place.  Only needs y and dy at t1 and t2.
Input
- state::AbstractState -- state of a stepper at time t2
- step_prev::Step -- solution at time t1 respectively
- tout -- time of requested output
- yout -- inplace y output
Ref: Hairer & Wanner p.190

TODO: tdir (I think this works for any order of t1 and t2 but needs
verifying.

TODO: fill dy

TODO: arbitrary order method (change step_prev::Step to step_prevs::Tuple{Step,N})

"""
function interpolate!{T,Y}(state::AbstractState,step_prev::Step{T,Y},tout::T,step_out::Step{T,Y})
    t1,y1,dy1 = output(step_prev)
    t2,y2,dy2 = output(state)
    if tout==t1
        copy!(step_out.y,y1)
    elseif tout==t2
        copy!(step_out.y,y2)
    else
        dt       = t2-t1
        theta    = (tout-t1)/dt
        for i=1:length(y1)
            step_out.y[i] =
                (1-theta)*y1[i] +
                theta*y2[i] +
                theta*(theta-1) *
                ( (1-2*theta)*(y2[i]-y1[i]) +
                  (theta-1)*dt*dy1[i] +
                  theta*dt*dy2[i])
        end
    end
    step_out.t = tout
    return nothing
end
