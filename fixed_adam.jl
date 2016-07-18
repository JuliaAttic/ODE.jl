#   Adapted from old AdamBashforth solver of ODE.jl:
#
#   ODE_MS  Solve non-stiff differential equations, multistep
#   fixed-step AdamBashforth method.
#
#   [T,X] = ODE_MS(ODEFUN, X0, TSPAN) with TSPAN = [T0:H:TFINAL]
#   integrates the system of differential equations x' = f(t,x) from time
#   T0 to TFINAL in steps of H with initial conditions X0. Function
#   ODEFUN(T,X) must return a column vector corresponding to f(t,x). Each
#   row in the solution array X corresponds to a time returned in the
#   column vector T.
#
#   function ode_ms(F, x0, tspan, order::Integer)
#       h = diff(tspan)
#       x = Array(typeof(x0), length(tspan))
#       x[1] = x0
#
#       if 1 <= order <= 4
#           b = ms_coefficients4
#       else
#           b = zeros(order, order)
#           b[1:4, 1:4] = ms_coefficients4
#           for s = 5:order
#               for j = 0:(s - 1)
#                   # Assign in correct order for multiplication below
#                   #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
#                   p_int = polyint(poly(diagm(-[0:j - 1; j + 1:s - 1])))
#                   b[s, s - j] = ((-1)^j / factorial(j)
#                                 / factorial(s - 1 - j) * polyval(p_int, 1))
#               end
#           end
#       end
#
#       # TODO: use a better data structure here (should be an order-element circ buffer)
#       xdot = similar(x)
#       for i = 1:length(tspan)-1
#           # Need to run the first several steps at reduced order
#           steporder = min(i, order)
#           xdot[i] = F(tspan[i], x[i])
#
#           x[i+1] = x[i]
#           for j = 1:steporder
#               x[i+1] += h[i]*b[steporder, j]*xdot[i-(steporder-1) + (j-1)]
#           end
#       end
#       return vcat(tspan), x
#   end

abstract FixedAdamSolver<:FixedStepSolver
immutable AdamBashforth{T}<:FixedAdamSolver
    tsteps::Vector{T}  # steps for the stepper to take
    order:: Int
    # TODO: make more general AbstractVector
end

immutable AdamMoulton{T}<:FixedAdamSolver
    tsteps::Vector{T}  # steps for the stepper to take
    order:: Int
    # TODO: make more general AbstractVector
end

function AdamBashforth{T,Y}(ode::IVP{T,Y}, opts::Dict)
    # :tsteps is required, all others are ignored...
    # TODO: figure out a way to make this constructor
    tsteps = similar(opts[:tsteps],T)
    copy!(tsteps, opts[:tsteps])

    if haskey(opts,"order")
      order = opts[:order]
    else
      order = 4
    end
    AdamBashforth{T}(tsteps,order)
end

function AdamMoulton{T,Y}(ode::IVP{T,Y}, opts::Dict)
    # :tsteps is required, all others are ignored...
    # TODO: figure out a way to make this constructor
    tsteps = similar(opts[:tsteps],T)
    copy!(tsteps, opts[:tsteps])

    if haskey(opts,"order")
      order = opts[:order]
    else
      order = 4
    end
    AdamMoulton{T}(tsteps,order)
end

# Note that the S-parameter serves no purpose here as this is a
# one-solver family (for now), thus it can be left off.  (But
# hard-code it in AState)
type AdamBashforthState{T,Y} <:AState{AdamBashforth, T, Y}
    dt::T     # (proposed) next time step

    t::T      # current time
    y::Y      # current solution
    dy::Y     # current derivative

    tpre::T  # time t-1
    ypre::Y  # solution at t-1
    dypre::Y  # solution at t-1
    prev_dys:: CircularBuffer{Y} # derivative at t-1

    step::Int # current step number
    finished::Bool # true if last step was taken

    b::Array{Rational{Int64},2}

end

type AdamMoultonState{T,Y} <:AState{AdamMoulton, T, Y}
    dt::T     # (proposed) next time step

    t::T      # current time
    y::Y      # current solution
    dy::Y     # current derivative

    tpre::T  # time t-1
    ypre::Y  # solution at t-1
    dypre::Y  # solution at t-1
    prev_dys:: CircularBuffer{Y} # derivative at t-1

    step::Int # current step number
    finished::Bool # true if last step was taken

    b_exp::Array{Rational{Int64},2}
    b_imp::Array{Rational{Int64},2}


end

#intialize the state
##Adam Bashforth Coefficients for Explicit Method
const ms_coefficients4 = Rational{Int64}[1        0        0         0
                                         -1//2    3//2     0         0
                                         5//12    -4//3    23//12    0
                                         -9//24   37//24   -59//24   55//24]
 ##Adam Moulton Coefficients for Implicit Method
 const am_imp_coefficients3 =Rational{Int64}[1         0        0         0
                                             1//2      1//2     0         0
                                             -1//12    8//12    5//12     0
                                             1//24     -5//24   19//24    9//24]


@compat function (::Type{AdamBashforthState}){M<:ExplicitODE,S<:AdamBashforth}(odep::Problem{M,S})
    @unpack odep: ode, solver
    y0 = ode.y0

    T = eltype(ode.tspan)
    Y = typeof(y0)
    EY = eltype(Y)

    tsteps = solver.tsteps
    @assert tsteps[1]==ode.tspan[1]
    @assert tsteps[end]==ode.tspan[2]
    order = solver.order

    N = length(tsteps)
    dof = length(y0)

    dt = tsteps[2]-tsteps[1]
    t = ode.tspan[1]

    y = deepcopy(y0)
    dy = similar(y0)
    ode.F!(t, y, dy)

    tpre = NaN
    ypre = zeros(y0)
    dypre = zeros(y0)
    prev_dys = CircularBuffer{Y}(order)
    push!(prev_dys, dy)

    step = 1
    finished = false

    # matrix which holds Adam Bashforth beta coefficients
    if 1 <= order <= 4
        b = ms_coefficients4
    else
        b = zeros(order, order)
        b[1:4, 1:4] = ms_coefficients4
        for s = 5:order
            for j = 0:(s - 1)
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
                p_int = polyint(poly(diagm(-[0:j - 1; j + 1:s - 1])))
                b[s, s - j] = ((-1)^j / factorial(j)
                               / factorial(s - 1 - j) * polyval(p_int, 1))
            end
        end
    end

    return AdamBashforthState(dt,
                   t, y, dy,
                   tpre, ypre, dypre, prev_dys,
                   step, finished,
                   b)
end

@compat function (::Type{AdamMoultonState}){M<:ExplicitODE,S<:AdamMoulton}(odep::Problem{M,S})
    @unpack odep: ode, solver
    y0 = ode.y0

    T = eltype(ode.tspan)
    Y = typeof(y0)
    EY = eltype(Y)

    tsteps = solver.tsteps
    @assert tsteps[1]==ode.tspan[1]
    @assert tsteps[end]==ode.tspan[2]
    order = solver.order

    N = length(tsteps)
    dof = length(y0)

    dt = tsteps[2]-tsteps[1]
    t = ode.tspan[1]

    y = deepcopy(y0)
    dy = similar(y0)
    ode.F!(t, y, dy)

    tpre = NaN
    ypre = zeros(y0)
    dypre = zeros(y0)
    prev_dys = CircularBuffer{Y}(order)
    push!(prev_dys, dy)

    step = 1
    finished = false

    # matrix which holds AdamBashforth and AdamMoulton beta coefficients
    if (1 <= order <= 4)
        b_imp = am_imp_coefficients3
        b_exp = ms_coefficients4
    else
        #calculating higher order coefficients for implicit Adam Moulton method
        b_imp = zeros(order, order)
        b_imp[1:4,1:4] = am_imp_coefficients3
        k = order - 1 # For explicit method, order = k+1
        for s = 4:k
            for j = 0:s
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(-1:s-1), except -(j-1), as roots)
                p_int = Polynomials.polyint(Polynomials.poly(diagm(-[-1:j - 2; j:s-1])))
                b_imp[s+1, s+1 - j] = ((-1)^j / factorial(j)
                               / factorial(s - j) * Polynomials.polyval(p_int, 1))
            end
        end
        b_exp = zeros(order, order)
        b_exp[1:4, 1:4] = ms_coefficients4
        k = order # For implicit method, order = k
        for s = 5:k
            for j = 0:(s - 1)
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
                p_int = Polynomials.polyint(Polynomials.poly(diagm(-[0:j - 1; j + 1:s - 1])))
                b_exp[s, s - j] = ((-1)^j / factorial(j)
                               / factorial(s - 1 - j) * Polynomials.polyval(p_int, 1))
            end
        end
    end

    return AdamMoultonState(dt,
                   t, y, dy,
                   tpre, ypre, dypre, prev_dys,
                   step, finished,
                   b_exp, b_imp)
end

# iteration
init{M,S<:AdamBashforth}(odep::Problem{M,S}) = AdamBashforthState(odep)
init{M,S<:AdamMoulton}(odep::Problem{M,S}) = AdamMoultonState(odep)

function onestep!(odep::ProblemFixedStep, st::AdamBashforthState)
    status = nothing # TODO
    F! = odep.ode.F!
    tsteps = odep.solver.tsteps
    order = odep.solver.order
    @unpack st: y, dy, t, dt, step, b, prev_dys
    steporder = min(step,order)
    dof = length(y)
    # integration for next y value
    ynext = y
    for j = 1:steporder
        # TODO: have dt be a circular buffer for dt for noneven fixed steps
        ynext += dt*b[steporder, j]*prev_dys[(step-(steporder-1) + (j-1)-1)%capacity(prev_dys)+1]
    end
    # update state
    ## previous data
    st.tpre = t
    st. ypre = y
    st.dypre = dy
    dynext = similar(y)
    F!(t+dt, ynext, dynext)
    push!(st.prev_dys,dynext)

    # current data
    st.t = t+dt
    st.dt = tsteps[st.step+1]-tsteps[st.step]
    st.step += 1
    st.y = ynext
    st.dy = dynext

    if st.step==length(odep.solver.tsteps)
        st.finished = true
    end

    return Status()
end


function onestep!(odep::ProblemFixedStep, st::AdamMoultonState)
    status = nothing # TODO
    F! = odep.ode.F!
    tsteps = odep.solver.tsteps
    order = odep.solver.order
    @unpack st: y, dy, t, dt, step, b_imp,b_exp, prev_dys
    steporder = min(step,order)
    dof = length(y)

    # integration for next y value
    ##(P)redict ynext using explicit Adam Bashforth coefficients
    ynext = y
    for j=1:steporder
        ynext += dt*b_exp[steporder, j]*prev_dys[(step-(steporder-1) + (j-1)-1)%capacity(prev_dys)+1]
    end

    ##(E)valuate function F at the approximate point t[i+1],ynext
    dynext = similar(y)
    F!(t+dt,ynext,dynext)
    push!(prev_dys,dynext)

    ##(C)orrect the formula using implicit Adam Moulton coefficients
    ynext = y
    for j = 1 : steporder
        ynext += dt*b_imp[steporder, j]*prev_dys[(step - (steporder -1) + j-1)%capacity(prev_dys)+1]
    end

    ##(E)valulate the function anew at corrected approximatiion t[i+1],ynext
    F!(t+dt, ynext,dynext)

    # update state
    ## previous data
    st.tpre = t
    st. ypre = y
    st.dypre = dy
    push!(st.prev_dys,dynext)

    # current data
    st.t = t+dt
    st.dt = tsteps[st.step+1]-tsteps[st.step]
    st.step += 1
    st.y = ynext
    st.dy = dynext

    if st.step==length(odep.solver.tsteps)
        st.finished = true
    end

    return Status()
end
