# ODE23S  Solve stiff systems based on a modified Rosenbrock triple
# (also used by MATLAB's ODE23s); see Sec. 4.1 in
#
# [SR97] L.F. Shampine and M.W. Reichelt: "The MATLAB ODE Suite," SIAM Journal on Scientific Computing, Vol. 18, 1997, pp. 1–22
#
# supports keywords: points = :all | :specified (using dense output)
#                    jacobian = G(t,y)::Function | nothing (FD)

# constants
const d = 1/(2 + sqrt(2))
const e32 = 6 + sqrt(2)

"""
Rosenbrock 23s solver, adapted from ode23s.

The different subtypes of Rosenbrock23s use different onestep! implementation.

- :v1 -- simple/naive port of ode23s, uses out-of-place functions

- :v2 -- advanced port: hooks into the trialstep!, rollback!,
         errorcontrol! framework.  Uses inplace function. TODO

"""
@with_kw immutable Rosenbrock23s{Name, N<:Number, Norm}<:AdaptiveSolver @deftype N
    # options
    reltol = 1.0e-5
    abstol = 1.0e-8
    points::Symbol=:all
    norm::Norm=Base.norm
    minstep=eps() # =abs(tspan[end] - tspan[1])/1e18
    maxstep=convert(N,Inf) # =abs(tspan[end] - tspan[1])/2.5
    initstep = 1e-8 # TODO
end
@compat function (::Type{Rosenbrock23s{:v1}}){T}(ode::IVP{T},
                                                   opts::Dict=Dict{Symbol,Any}())
    _construct_Rosenbrock23s(T, Rosenbrock23s{:v1}, ode, opts)
end
@compat function (::Type{Rosenbrock23s{:v2}}){T}(ode::IVP{T},
                                                   opts::Dict=Dict{Symbol,Any}())
    _construct_Rosenbrock23s(T, Rosenbrock23s{:v2}, ode, opts)
end
function _construct_Rosenbrock23s(T, RB23s, ode, opts)
    Norm = typeof(Base.norm)
    rb = RB23s{T, Norm}()
    args = []
    for f in fieldnames(Rosenbrock23s)
        if haskey(opts, f)
            push!(args, opts[f])
        else
            push!(args, getfield(rb, f))
        end
        if f==:norm
            Norm = typeof(args[end])
        end
    end
    RB23s{T, Norm}(args...)
end

type Rosenbrock23sState{T,Y,FF,JJ,Y2,LU} <: AState{Rosenbrock23s,T,Y}
    dt::T     # (proposed) next time step size
    t::T      # current time
    y::Y      # current solution
    dy::Y     # current derivative

    tpre::T  # time t-1
    ypre::Y  # solution at t-1
    dypre::Y # derivative at t-1

    step::Int # current step number
    finished::Bool # true when done with last step
    staleJ::Bool # whether to recalculate jacobian

    # hacks
# ob: what are these used for ?
    F::FF # out-of-place function
    J::JJ # out-of-place function

    # work arrays
    k1::Y
    k2::Y
    k3::Y
    jac::Y2 # Jacobian (and lufact! storage)
    W::LU  # and its lufact

end
function Rosenbrock23sState{M<:ExplicitODE,S<:Rosenbrock23s}(odep::Problem{M,S})
    @unpack odep: ode, solver
    @unpack ode: tspan, y0, F!, G!, J!
    @assert (G!)==nothing "Cannot use implicit ODEs or mass-matrix with this solver"

    F = outofplace(F!)

    T = eltype(tspan)
    Y = typeof(y0)
    EY = eltype(Y)

    dt::T = solver.initstep
    t = tspan[1]
    y = deepcopy(y0)
    dy = similar(y0)
    ode.F!(t, y, dy)

    tpre::T = NaN
    ypre = zeros(y0)
    dypre = zeros(y0)

    step = 1
    finished = false

    # deal with Jacobian
    if (J!)==nothing
        # fallback finite-difference
        J = fd_jacobian(F)
    else
        J = outofplace(J!)
    end
    jac = J(t,y)
    jac = eye(jac) - dt*d*jac    # get Jacobian of F wrt y
    W = lufact!(jac)
    staleJ = false

    return Rosenbrock23sState(dt, t, y, dy,
                              tpre, ypre, dypre,
                              step, finished, staleJ,
                              F,J,
                              zeros(y0),zeros(y0),zeros(y0),
                              jac, W
                              )
end

init{M,S<:Rosenbrock23s}(odep::Problem{M,S}) = Rosenbrock23sState(odep)

#####
# Naive implementation :v1
#####

# Rosenbrock23s{:v1} onestep! implementation
function onestep!{M,S<:Rosenbrock23s{:v1}}(odep::ProblemAdaptiveStep{M,S}, st::Rosenbrock23sState)
    @unpack st: dt, t, y, dy, ypre, dypre, F, J
    @unpack odep: ode
    @unpack ode: tspan
    @unpack odep.solver: maxstep, minstep, abstol, reltol

    # use as work-arrays:
    ynext = ypre
    dynext = dypre


    jac = J(t,y)    # get Jacobian of F wrt y
    # note: if there is a mass matrix M on the lhs of the ODE, i.e.,
    #   M * dy/dt = F(t,y)
    # we can simply replace eye(jac) by M in the following expression
    # (see Sec. 5 in [SR97])
    W = lufact( eye(jac) - dt*d*jac )

    accept = false
    while !accept

        F0 = dy

        # approximate time-derivative of F
        T = dt*d*(F(t + dt/100, y) - F0)/(dt/100)

        # modified Rosenbrock formula
        k1 = W\(F0 + T)
        F1 = F(t + 0.5*dt, y + 0.5*dt*k1)
        k2 = W\(F1 - k1) + k1
        ynext = y + dt*k2
        F2 = F(t + dt, ynext)
        k3 = W\(F2 - e32*(k2 - F1) - 2*(k1 - F0) + T )

        err = (abs(dt)/6)*norm(k1 - 2*k2 + k3) # error estimate
        delta = max(reltol*max(norm(y),norm(ynext)), abstol) # allowable error

        # check if ynext solution is acceptable
        if  err <= delta
            accept = true
            dynext = F2
            t += dt
        end

        # update of the step size
        dt = tdir(ode)*min( maxstep, abs(dt)*0.8*(delta/err)^(1/3) )

        # don't overstep tspan[2]:
        if tdir(ode)*(tspan[2]-t) <= tdir(ode)*dt
            dt = tspan[2]-t
        end
        if dt!=0 && abs(dt)<max(minstep,eps(t)) # TODO: use this generally
            error("Minimum step size $(max(minstep,eps(t))) reached")
        end
    end
    st.finished = dt==0

    # update state
    st.dt = dt

    st.tpre = st.t
    st.ypre = y
    st.dypre = dy

    st.t = t
    st.step += 1

    st.y = ynext
    st.dy = dynext

    return Status()
end


#####
# Better implementation :v2
#####

# This takes one step of size dt.
#
# TODO: make use of in-place functions, at least where possible.
function trialstep!{M,S<:Rosenbrock23s{:v2}}(odep::ProblemAdaptiveStep{M,S}, st::Rosenbrock23sState)
    @unpack st: dt, t, y, dy, ypre, dypre, F, J, jac, W, staleJ

    # Use the *pre as work-arrays and leave y alone until the state gets updated.
    ynext = ypre
    dynext = dypre

    if staleJ
        jac[:] = eye(jac) - dt*d*J(t,y)    # get Jacobian of F wrt y
        # Note: if there is a mass matrix M on the lhs of the ODE, i.e.,
        #   M * dy/dt = F(t,y)
        # we can simply replace eye(jac) by M in the following expression
        # (see Sec. 5 in [SR97])
        st.W = lufact!(jac)
        W = st.W
    end

    F0 = dy

    # approximate time-derivative of F
    # TODO: use ForwardDiff?
# ob: consolating ForwardDiff documentation http://www.juliadiff.org/ForwardDiff.jl/api.html#derivatives-of
# dFdt = ForwardDiff.derivative(F,t)
    dFdt = dt*d*(F(t + dt/100, y) - F0)/(dt/100)

    # modified Rosenbrock formula
    k1 = W\(F0 + dFdt)
    F1 = F(t + 0.5*dt, y + 0.5*dt*k1)
    k2 = W\(F1 - k1) + k1

    ynext = y + dt*k2
    dynext = F(t + dt, ynext)
    F2 = dynext
    k3 = W\(F2 - e32*(k2 - F1) - 2*(k1 - F0) + dFdt )

    # update state
    st.k1 = k1
    st.k2 = k2
    st.k3 = k3

    st.tpre = t
    st.t = t+dt
    st.ypre = y
    st.dypre = dy
    st.y = ynext
    st.dy = dynext
    st.staleJ = true
    return nothing # some stats
end

# Calculates the error and new dt.  Sets the finished flag if the last
# step was successful and up to tspan[2].
function errorcontrol!{M,S<:Rosenbrock23s{:v2}}(odep::ProblemAdaptiveStep{M,S}, st::Rosenbrock23sState)
    # TODO: make a generic error control method
    @unpack st: k1, k2, k3, dt, t, tpre, y, ypre
    @unpack odep.solver: abstol, reltol, maxstep, minstep, norm
    @unpack odep: ode
    @unpack odep.ode: tspan
    err = (abs(dt)/6)*norm(k1 - 2*k2 + k3) # error estimate
    delta = max(reltol*max(norm(y),norm(ypre)), abstol) # allowable error

    err_out = err/delta
    if isnan(err_out)
        error("NaN error encountered")
    end

    dt = tdir(ode)*min( maxstep, abs(dt)*0.8*(delta/err)^(1/3) )
    if err_out>=1 # step is rejected and t is rewound to tpre
        if tdir(ode)*(tspan[2]-tpre) <= tdir(ode)*dt
            dt = tspan[2]-tpre
        end
        if dt==0
            error("This shouldn't happen")
        end
    else # step is accepted
        if tdir(ode)*(tspan[2]-t) <= tdir(ode)*dt
            dt = tspan[2]-t
        end
        st.finished = dt==0
    end

    # update state
    st.dt = dt
    return err_out, nothing
end

# If the step was rejected, this rolls-back the state variable such
# that another trial step can be taken from the old t.

# ob: I guess rollback!() isn't integrated yet? Is the idea behind rollback!()
#     that retrying a step may not require as much recomputation is just redoing
#     the step with a smaller dt?

function rollback!{M,S<:Rosenbrock23s{:v2}}(odep::ProblemAdaptiveStep{M,S}, st::Rosenbrock23sState)
    # just load the values from the *pre back
    st.y = st.ypre
    st.dy = st.dypre
    st.t = st.tpre
    st.finished = false
    st.staleJ = false
    return nothing
end
