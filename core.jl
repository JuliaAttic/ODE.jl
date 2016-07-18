# This provides the core of ODE.jl into which the solvers hook into
# from one side and the user interfaces from the other.  It is based
# around iteration and much inspired by @pwl's PR 49.

#####
# ODE/DAE problems: mathematical formulation
#####
"""
Defines the mathematical part of a ODE/DAE specified in the general
form:

`F(t, y) =  G(t, y, dy) with y(t0)= y0`

although not all solvers will support any combinations of `F` and `G`.
Note that not specifying `G` amounts to `G=dy/dt`.


`T<:Number` is the type of the time, `Y<:AbstractVector` is the type
of `y`. `F`,`G` and `J` are the types of the functions or the latter
two can also be matrices.

- `tspan` -- tuple `(start_t,end_t)`
- `y0` -- initial condition
- `F!` -- in-place `F` function `F!(t,y,res)`.  If `F=0` set to `nothing`.
- `G!` -- in-place `G` function `G!(t,y,dy,res)`.  If `G=dy/dt` then
          set it to `nothing` (or `dy` if the solver supports this).  Can
          also be a mass matrix for a RHS `M dy/dt`
- `J!` -- in-place Jacobian function `J!(t,y,dy,res)`.  Or a sparsity
          pattern of Jacobian used (in the future) for
          finite-difference or auto-diff Jacobians.  (c.f. matrix
          coloring, TODO)
"""
immutable IVP{T<:Number, Y<:AbstractVector, F, G, J}
    tspan  ::Tuple{T,T}
    y0  ::Y
    # Note, in Julia 0.5, each function needs its own type for type
    # stable higher order functions taking IVP as argument.
    F!  ::F
    G!  ::G # can also contain a mass-matrix
    J!  ::J # can also contain a sparsity matrix
end

# Using typealias, different types of ODE/DAEs can be encoded:

typealias ExplicitODE{T,Y,F,J}   IVP{T,Y,F,Void,J}
@compat (::Type{ExplicitODE}){T,Y,F}(tspan::Tuple{T,T}, y0::Y, F!::F) = IVP{T,Y,F,Void,Void}(tspan,y0,F!,nothing,nothing)

typealias ExplicitODEwithMass{T,Y,F,M<:AbstractMatrix,J}   IVP{T,Y,F,M,J}
@compat (::Type{ExplicitODEwithMass}){T,Y,F,M}(tspan::Tuple{T,T}, y0::Y, F!::F, mass::M) = IVP{T,Y,F,M,Void}(tspan,y0,F!,mass,nothing)

typealias ImplicitODE{T,Y,G,J}   IVP{T,Y,Void,G,J}
@compat (::Type{ImplicitODE}){T,Y,G,J}(tspan::Tuple{T,T}, y0::Y, G!::G, J!::J) = IVP{T,Y,Void,G,J}(tspan,y0,nothing,G!,J!)

typealias IMEX_ODE{T,Y,F,G,J}   IVP{T,Y,F,G,J}
# no constructor needed

"Return end time of integration"
finaltime(ode::IVP) = ode.tspan[end]

"Type of time"
t_type{T}(::Type{IVP{T}}) = T
t_type{T}(::IVP{T}) = T

"Type of y (will be <:AbstractVector)"
y_type{T,Y}(::Type{IVP{T,Y}}) = Y
y_type{T,Y}(::IVP{T,Y}) = Y

#####
# Solvers
#####

"""
Numerical ODE/DAE solver.

The subtypes hold the options for the solver.  Subtype either
AdaptiveSolver or FixedStepSolver.

A solver instance also holds all the options for a problem.  This
means that the type parameters of a subtype of ASolver relate to the
types of the solver options.  Thus, in general, a IVP needs to be
supplied to instantiate the solver, as the type parameters are
dependent on the IVP parameters.  Therefore, define a constructor:

`ASolverSubT(ode::IVP, opts::Dict{Symbol}=Dict{Symbol,Any}()) = ...`
"""
abstract ASolver
name{S<:ASolver}(s::Type{S}) = s.name.name
name(s::ASolver) = name(typeof(s))
abstract AdaptiveSolver <: ASolver
abstract FixedStepSolver <: ASolver

"""
Returns a dictionary of supported options for a specific IVP &
Type{Solver} pair.
"""
get_solver_defopts{S<:ASolver}(ode::IVP, Solver::Type{S}) = type2dict(Solver(ode))

order(s::ASolver) = error("This method needs implementing by solver $s.")

"""
Holds the solver status after onestep.

TODO
"""
immutable Status end


#####
# Full problem consists of IVP + ASolver
#####

"""
Problem holds a ODE/DAE problem, which consists of

- mathematical definition `IVP`
- numerical solver `ASolver`, including its options
- dictionary of all user options (say also for dense output, events, etc.)

TODO: it would be cool to have a diagnostic tool to check which
options end up being used.  Both to catch typos and to just find out.
"""
immutable Problem{EQ<:IVP, S<:ASolver}
    ode::EQ
    solver::S
    opts::Dict{Symbol,Any}
end
function Problem{EQ}(ode::EQ, Solver::DataType, opts::Dict{Symbol}=Dict{Symbol,Any}())
    solver = Solver(ode, opts)
    Problem{EQ,typeof(solver)}(ode,solver,opts)
end

typealias ProblemFixedStep{EQ,S<:FixedStepSolver} Problem{EQ,S}
typealias ProblemAdaptiveStep{EQ,S<:AdaptiveSolver} Problem{EQ,S}

Base.length(odep::ProblemFixedStep) = length(odep.solver.tsteps)

"Return end time of integration"
finaltime(odep::Problem) = finaltime(odep.ode)

# Type of time
t_type{EQ}(::Type{Problem{EQ}}) = t_type(EQ)
t_type{EQ}(::Problem{EQ}) = t_type(EQ)

# Type of y (will be <:AbstractVector)
y_type{EQ}(::Type{Problem{EQ}}) = y_type(EQ)
y_type{EQ}(::Problem{EQ}) = y_type(EQ)


#####
# Iteration of a `Problem`
#####

"""
Subtypes of `AState` are used to hold the state of a solver.  Given
this state it needs to be possible to:

- take the next step
- create dense output between the current time and the last time

Parameters are:

- S -- a solver.
- T -- type of time.  Probably `<:AbstractFloat`
- Y -- type of y0.  `<:AbstractVector` (note that scalar like types need
       to be wrapped in a 1-element vector)
"""
abstract AState{S, T, Y}
# was:
# abstract AState{S<:ASolver, T<:Number, Y<:AbstractVector}
# Note, the constraints need repeating in subtypes!
# https://github.com/JuliaLang/julia/issues/9441
# thus leave them off for now as it leads to tricky errors

# Type of time
t_type{S,T,Y}(::Type{AState{S,T,Y}}) = T
t_type{AS<:AState}(::Type{AS}) = t_type(supertype(AS))
t_type(st::AState) = t_type(supertype(typeof(st)))

# Type of y (will be <:AbstractVector)
y_type{S,T,Y}(::Type{AState{S,T,Y}}) = Y
y_type{AS<:AState}(::Type{AS}) = y_type(supertype(AS))
y_type(st::AState) = y_type(supertype(typeof(st)))

"Return solver belonging to a state variable/type"
getsolver{S}(::Type{AState{S}}) = S
getsolver{AS<:AState}(T::Type{AS}) = getsolver(supertype(T))
getsolver{S}(::AState{S}) = S

# Interface to get to fields:
gety(st::AState) = st.y
getdy(st::AState) = st.dy
gett(st::AState) = st.t
gety_pre(st::AState) = st.ypre
getdy_pre(st::AState) = st.dypre
gett_pre(st::AState) = st.tpre

# TODO an interface like so might be nice, especially for multistep
# solvers.  Could use a circular buffer:
# gety(st::AState,i) = st.ys[i]
# getdy(st::AState,i) = st.dys[i]
# gett(st::AState,i) = st.ts[i]

# Iteration: take one step on a ODE/DAE `Problem`
#
# Define:
# start(iter) -> state
# next(iter, state) -> item, state
# done(iter, state) -> bool
#
# Note that `item` and `state` are the same.
#
# Note defining `length` maybe a bad idea as collect does not work...
Base.start(odep::Problem) = init(odep)
Base.next(odep::Problem, st::AState) = (status = onestep!(odep, st); ((status,st),st))
Base.done(odep::Problem, st::AState) = finished(st)

#####
# Interface to implement by solvers to hook into iteration
#####
#
# See runge_kutta.jl and rosenbrock.jl for example implementations.

"Tests for finished."
finished(st::AState) = st.finished

"""
Take one step.  This is the core function to be implemented by a
solver.  Note that adaptive solvers may want to implement only some of
the substeps.
"""
onestep!(odep::ProblemFixedStep, st::AState) =
    error("Function `onestep!` needs to be implemented for fixed step solver $S")

# This onestep! implementation breaks the looping down into different
# parts.  This currently works with the second Rosenbrock23s solver.
#
# TODO: This be used with fixed-step solvers with errorcontrol!()==0.
# Would this be clearer?
function onestep!(odep::ProblemAdaptiveStep, st::AState)
    accept = false
    while !accept
        stats = trialstep!(odep, st)
        err, stats = errorcontrol!(odep, st)
        if err<=1
            accept = true
        else
            rollback!(odep, st)
        end
    end
    return Status()
end

# ob: In response to your above question, I think a simple comment that explains
# how fixed step functions have no need for trialsteps would be clear enough.

"""
Advances the solution to new state by a given time step.  Updates
state in-place such that it reflects the new state.

Returns the stats for this step (TODO).
"""
trialstep!{EQ,S<:AdaptiveSolver}(odep::Problem{EQ,S}, st::AState) =
    error("Function `trialstep!` and companions (or alternatively `onestep!`) need to be implemented for adaptive solver $S")

"""
Reverts (in-place) the state back to the previous state after a failed
trial step.  The reverting needn't be 100% as long as a new trial step
can be calculated from it.

Returns nothing.
"""
rollback!{EQ,S<:AdaptiveSolver}(odep::Problem{EQ,S}, st::AState) =
    error("Function `rollback!` and companions (or alternatively `onestep!`) need to be implemented for adaptive solver $S")


"""
Estimates the error (such that a step is accepted if err<=1), a new dt
and a new order.  Updates state with new dt and order (as appropriate).

Returns err & stats (TODO).
"""
errorcontrol!{EQ,S<:AdaptiveSolver}(odep::Problem{EQ,S}, st::AState) =
    error("Function `errorcontrol!` and companions (or alternatively `onestep!`) need to be implemented for adaptive solver $S")
