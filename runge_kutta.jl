# Adapted from old RK4 solver of ODE.jl:
#
#ODE4  Solve non-stiff differential equations, fourth order
#   fixed-step Runge-Kutta method.
#
#   [T,X] = ODE4(ODEFUN, X0, TSPAN) with TSPAN = [T0:H:TFINAL]
#   integrates the system of differential equations x' = f(t,x) from time
#   T0 to TFINAL in steps of H with initial conditions X0. Function
#   ODEFUN(T,X) must return a column vector corresponding to f(t,x). Each
#   row in the solution array X corresponds to a time returned in the
#   column vector T.
# function ode4(F, x0, tspan)
#     h = diff(tspan)
#     x = Array(typeof(x0), length(tspan))
#     x[1] = x0
#
#     midxdot = Array(typeof(x0), 4)
#     for i = 1:length(tspan)-1
#         # Compute midstep derivatives
#         midxdot[1] = F(tspan[i],         x[i])
#         midxdot[2] = 2*F(tspan[i]+h[i]./2, x[i] + midxdot[1].*h[i]./2)
#         midxdot[3] = 2*F(tspan[i]+h[i]./2, x[i] + midxdot[2].*h[i]./2)
#         midxdot[4] = F(tspan[i]+h[i],    x[i] + midxdot[3].*h[i])
#
#         # Integrate
#         x[i+1] = x[i] + 1/6 .*h[i].*sum(midxdot)
#     end
#     return [tspan], x
# end

abstract ExplicitRK<:FixedStepSolver
immutable RK4{T}<:ExplicitRK
    tsteps::Vector{T}  # steps for the stepper to take
    # TODO: make more general AbstractVector
end
function RK4{T,Y}(ode::IVP{T,Y}, opts::Dict)
    # :tsteps is required, all others are ignored...
    # TODO: figure out a way to make this constructor
    tsteps = similar(opts[:tsteps],T)
    copy!(tsteps, opts[:tsteps])
    RK4{T}(tsteps)
end

# Note that the S-parameter serves no purpose here as this is a
# one-solver family (for now), thus it can be left off.  (But
# hard-code it in AState)
type RKState{T,Y} <:AState{RK4, T, Y}
    dt::T     # (proposed) next time step

    t::T      # current time
    y::Y      # current solution
    dy::Y     # current derivative

    tpre::T  # time t-1
    ypre::Y  # solution at t-1
    dypre::Y # derivative at t-1

    step::Int # current step number
    finished::Bool # true if last step was taken

    # work arrays
    mid_y::Y
    mid_dy::Vector{Y}
end
function RKState{M<:ExplicitODE,S<:RK4}(odep::Problem{M,S})
    @unpack odep: ode, solver
    y0 = ode.y0

    T = eltype(ode.tspan)
    Y = typeof(y0)
    EY = eltype(Y)

    tsteps = solver.tsteps
    @assert tsteps[1]==ode.tspan[1]
    @assert tsteps[end]==ode.tspan[2]

    N = length(tsteps)
    dof = length(y0)

    dt = tsteps[2]-tsteps[1]
    t = ode.tspan[1]

    y = deepcopy(y0)
    dy = similar(y0)
    ode.F!(t, y, dy)

    # ob: I thought including these assignments might make it a little clearer
    # to a user what is going on. Though, it is longer.
    tpre = NaN
    ypre = zeros(y0)
    dypre = zeros(y0)

    step = 1
    finished = false

    mid_y = similar(y0)
    mid_dy = Y[similar(y0), similar(y0), similar(y0)]
    return RKState(dt,
                   t, y, dy,
                   tpre, ypre, dypre,
                   step, finished,
                   mid_y, mid_dy)
end

const onestep_counter_RK4 = [0]
# iteration
init{M,S<:RK4}(odep::Problem{M,S}) = RKState(odep)
function onestep!(odep::ProblemFixedStep, st::RKState)
    onestep_counter_RK4[1] += 1
    status = nothing # TODO
    F! = odep.ode.F!
    tsteps = odep.solver.tsteps
    @unpack st: y, dy, t, dt, mid_y, mid_dy

    ynext = st.ypre
    dynext = st.dypre
    dof = length(y)

    for i in 1:dof
        mid_y[i] = y[i] + dy[i]*dt/2
    end
    F!(t+dt/2, mid_y, mid_dy[1])

    for i in 1:dof
        mid_y[i] = y[i] + mid_dy[1][i]*dt
    end
    F!(t+dt/2, mid_y, mid_dy[2])

    for i in 1:dof
        mid_y[i] = y[i] + mid_dy[2][i]*dt
    end
    F!(t+dt  , mid_y, mid_dy[3])

    # Integrate
    for i in eachindex(st.ypre)
        ynext[i] = y[i] + 1/6 .*dt.*(dy[i] +
                                    2*mid_dy[1][i] +
                                    2*mid_dy[2][i] +
                                    mid_dy[3][i])
    end

    # update state
    F!(t+dt, ynext, dynext)


    st.t = t+dt
    st.dt = tsteps[st.step+1]-tsteps[st.step]
    st.step += 1

    st.tpre = t
    st.ypre = y
    st.dypre = dy

    st.y = ynext
    st.dy = dynext

    if st.step==length(odep.solver.tsteps)
        st.finished = true
    end

    return Status()
end
