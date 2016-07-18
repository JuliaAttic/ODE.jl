# This file contains the implementation of explicit Runkge-Kutta
# solver from (Hairer & Wanner 1992 p.134, p.165-169).
#
# Adapted from PR49


"""
An adaptive Runge-Kutta solver

TODO
"""
@with_kw immutable RKAdaptive{Name, N<:Number, Norm}<:AdaptiveSolver @deftype N
    tableau::Tableau{}
    # options
    reltol = 1.0e-5
    abstol = 1.0e-8
    points::Symbol=:all
    norm::Norm=Base.norm
    minstep=eps() # =abs(tspan[end] - tspan[1])/1e18
    maxstep=convert(N,Inf) # =abs(tspan[end] - tspan[1])/2.5
    initstep = 1e-8 # TODO
end
"""
An fixed step Runge-Kutta solver
"""
immutable RKFixedStep{Name, T}<:FixedStepSolver
    tableau::Tableau{T}
    tsteps::Vector{T}
end
@compat function (::Type{RKFixedStep{Name}}){Name,T}(ode::IVP{T}, opts::Dict)
    # :tsteps is required, all others are ignored...
    tsteps = similar(opts[:tsteps],T)
    copy!(tsteps, opts[:tsteps])

    tab = convert(TableauRKExplicit{T} ,tableaus_rk_explicit[Name])
    RKFixedStep{Name, T}(tab,tsteps)
end

typealias RKSolver Union{RKAdaptive, RKFixedStep}

order(solver::RKSolver) = order(solver.tableau)[1]

####
# all solvers
rk_fixedstep_solvers = Dict{Symbol,DataType}()
for name_ in keys(tableaus_rk_explicit)
    rk_fixedstep_solvers[name_] = RKFixedStep{name_}
end

type RKState_2{S,T,Y} <: AState{S,T,Y}
    dt::T     # (proposed) next time step size
    t::T      # current time
    y::Y      # current solution
    dy::Y     # current derivative

    tpre::T  # time t-1
    ypre::Y  # solution at t-1
    dypre::Y # derivative at t-1

    step::Int # current step number
    finished::Bool # true when done with last step
    timeout::Int # timeout before step size increase is allowed

    # work arrays
    yerr::Y # in-place error
    ytmp::Y # used in calc_next_k!
    ks  ::Vector{Y}
end
function RKState_2{M<:ExplicitODE, S<:RKSolver}(odep::Problem{M,S})
    @unpack odep: ode, solver
    @unpack ode: y0, F!, tspan
    @unpack solver: tableau

    T = eltype(tspan)
    Y = typeof(y0)
    EY = eltype(Y)

    tsteps = solver.tsteps
    @assert tsteps[1]==tspan[1]
    @assert tsteps[end]==tspan[2]

    N = length(tsteps)
    dof = length(y0)

    if S<:AdaptiveSolver
        error("TODO: branch between fixed step and adaptive")
    elseif S<:FixedStepSolver
        nothing
    else
        error("Need AdaptiveSolver or FixedStepSolver")
    end

    dt = tsteps[2]-tsteps[1]
    t = tspan[1]

    y = deepcopy(y0)

    # ob: I thought including these assignments might make it a little clearer
    # to a user what is going on. Though, it is longer.
    tpre = NaN
    ypre = zeros(y0)
    ytmp = similar(y0)

    step = 1
    finished = false
    timeout = 0 # for step control

    # work
    yerr = similar(y)
    lk = lengthks(tableau)
    ks = Array(typeof(y0), lk) # ks
    for i = 1:lk # we have to allocate each component separately
        ks[i]=zero(y0)
    end
    # pre-initialize ks[1]
    F!(t,y,ks[1])

    # dy is in ks:
    dy = ks[1]
    dypre = ks[1] # just for now

    return RKState_2{S,T,Y}(dt, t, y, dy,
                      tpre, ypre, dypre,
                      step, finished, timeout,
                      yerr, ytmp, ks)
end
init{M,S<:RKSolver}(odep::Problem{M,S}) = RKState_2(odep)

#####################
# Fixed step method #
#####################
const onestep_counter_rk = [0]
function onestep!(odep::ProblemFixedStep, st::RKState_2)
    onestep_counter_rk[1] += 1
    @unpack st: step, ks, y, t, dt, ypre, ytmp
    @unpack odep.ode: F!
    @unpack odep.solver: tableau, tsteps
    @unpack tableau: b, isFSAL

    ynext = st.ypre
    dof  = length(y)

    # initialize
    copy!(ynext, y)

    for s=1:size(b,2)
        calc_next_k!(ks, ytmp, y, s, F!, t, dt, dof, tableau)
        for d=1:dof
            ynext[d] += dt * b[1,s]*ks[s][d]
        end
    end
    # update state
    st.y = ynext
    st.ypre = y
    st.t += dt

    st.dt = tsteps[st.step+1]-tsteps[st.step]

    st.t = t+dt
    st.dt = tsteps[st.step+1]-tsteps[st.step]
    st.step += 1

    st.tpre = t
    st.dypre = ks[1]

    st.y = ynext

    if isFSAL
        ks[1],ks[end] = ks[end],ks[1]
    else
        F!(st.t,st.y,ks[1])
    end

    if st.step==length(odep.solver.tsteps)
        st.finished = true
    end

    return Status()
end

function calc_next_k!{Ty}(ks::Vector, ytmp::Ty, y, s, fn, t, dt, dof, tableau)
    # Calculates the next ks and puts it into ks[s]
    # - ks and ytmp are modified inside this function.

    # Needed interface:
    # On components: +, *
    # On y0 container: setindex!, getindex, fn

    @unpack tableau: a,c

    copy!(ytmp, y)
    for ss=1:s-1, d=1:dof
        ytmp[d] += dt * ks[ss][d] * a[s,ss]
    end
    fn(t + c[s]*dt, ytmp, ks[s])
    nothing
end


########################
# Adaptive step method #
########################


# function next{RKSA<:RKStepperAdaptive}(sol::Solver{RKSA}, state::RKState)

#     const timeout_const = 5

#     # the initial values
#     dt      = state.dt          # dt is the previous stepisze, it is
#     # modified inside the loop
#     timeout = state.timeout
#     work    = state.work
#     step    = state.step
#     tableau = sol.stepper.tableau

#     # The while loop continues until we either find a stepsize which
#     # leads to a small enough error or the stepsize reaches
#     # prob.minstep

#     # trim the inital stepsize to avoid overshooting
#     dt      = min(dt, sol.options.tstop-state.step.t)

#     while true

#         # Do one step (assumes ks[1]==f0).  After calling work.ynew
#         # holds the new step.
#         # TODO: return ynew instead of passing it as work.ynew?

#         # work.y and work.yerr and work.ks are updated after this step
#         rk_embedded_step!(work, sol.ode, tableau, step, dt)

#         # changes work.yerr
#         err, newdt, timeout = stepsize_hw92!(work, step, tableau, dt, timeout, sol.options)

#         # trim again in case newdt > dt
#         newdt = min(newdt, sol.options.tstop-state.step.t)

#         if abs(newdt) < sol.options.minstep  # minimum step size reached, break
#             # passing the newdt to state will result in done()
#             state.dt = newdt
#             break
#         end

#         if err > 1 # error is too large, repeat the step with smaller dt
#             # redo step with smaller dt and reset the timeout
#             dt      = newdt
#             timeout = timeout_const
#         else
#             # step is accepted

#             # preload ks[1] for the next step
#             if sol.stepper.tableau.isFSAL
#                 copy!(work.ks[1],work.ks[end])
#             else
#                 sol.ode.F!(step.t+dt, work.ynew, work.ks[1])
#             end

#             # Swap bindings of y and ytrial, avoids one copy
#             step.y, work.ynew = work.ynew, step.y

#             # Update state with the data from the step we have just
#             # made:
#             step.t += dt
#             state.dt = newdt
#             state.timeout = timeout
#             break
#         end
#     end
#     return ((step.t,step.y),state)
# end


# ##########################
# # Lower level algorithms #
# ##########################

# function rk_embedded_step!(work      ::RKWorkArrays,
#                            ode       ::ExplicitODE,
#                            tableau   ::Tableau,
#                            last_step ::Step,
#                            dt)
#     # Does one embedded R-K step updating work.ynew, work.yerr and work.ks.
#     # Assumes that work.ks[:,1] is already calculated!
#     # Modifies work.y, work.ynew and work.yerr only

#     y      = last_step.y
#     dof    = length(y)
#     b      = tableau.b

#     fill!(work.ynew, zero(eltype(y)))
#     fill!(work.yerr, zero(eltype(y)))

#     for s=1:lengthks(tableau)
#         # we skip the first step beacause we assume that work.ks[1] is
#         # already computed
#         if s > 1
#             calc_next_k!(work, s, ode, tableau, last_step, dt)
#         end
#         for d=1:dof
#             work.ynew[d] += b[1,s]*work.ks[s][d]
#             work.yerr[d] += b[2,s]*work.ks[s][d]
#         end
#     end

#     for d=1:dof
#         work.yerr[d] = dt*(work.ynew[d]-work.yerr[d])
#         work.ynew[d] = y[d] + dt*work.ynew[d]
#     end

# end


# function stepsize_hw92!{T}(work,
#                            last_step ::Step,
#                            tableau   ::TableauRKExplicit,
#                            dt        ::T,
#                            timeout,
#                            options   ::Options)
#     # Estimates the error and a new step size following Hairer &
#     # Wanner 1992, p167 (with some modifications)
#     #
#     # If timeout>0 no step size increase is allowed, timeout is
#     # decremented in here.
#     #
#     # Returns the error, newdt and the number of timeout-steps
#     #
#     # TODO:
#     # - allow component-wise reltol and abstol?
#     # - allow other norms

#     ord = minimum(order(tableau))
#     timout_after_nan = 5
#     # fac = T[0.8, 0.9, (0.25)^(1/(ord+1)), (0.38)^(1/(ord+1))][1]
#     fac = T(8//10)
#     facmax = T(5) # maximal step size increase. 1.5-5
#     facmin = 1./facmax  # maximal step size decrease. ?
#     dof = length(last_step.y)

#     if findfirst(options.isoutofdomain,work.y) != 0
#         return T(10), dt*facmin, timout_after_nan
#     end

#     # in-place calculate yerr./tol
#     for d=1:dof
#         y0 = last_step.y[d] # TODO: is this supposed to be the last successful step?
#         y1 = work.ynew[d]    # the approximation to the next step
#         sci = (options.abstol + options.reltol*max(norm(y0),norm(y1)))
#         work.yerr[d] ./= sci # Eq 4.10
#     end

#     # TOOD: should we use options.norm here as well?
#     err   = norm(work.yerr) # Eq. 4.11
#     newdt = min(options.maxstep, dt*max(facmin, fac*(1/err)^(1/(ord+1)))) # Eq 4.13 modified

#     if timeout > 0
#         newdt = min(newdt, dt)
#         timeout -= 1
#     end

#     return err, newdt, timeout
# end


# # Calculates the next ks and puts it into ks[s]
# # - ks and ytmp are modified inside this function.
# function calc_next_k!(ks, ytmp, y, i, F!, t, dt, dof, btab)
#     @unpack tableau: a,c
#     dof = length(y)
#     t, a, c = last_step.t, tableau.a, tableau.c

#     copy!(y,last_step.y)
#     for j=1:i-1
#         for d=1:dof
#             y[d] += dt * ks[j][d] * a[i,j]
#         end
#     end
#     ode.F!(t + c[i]*dt, y, ks[i])
# end
