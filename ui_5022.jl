# High(er)-level user interface
"""
Converts the inputted y0,tsteps,tout into good types, if needed.

TODO: implement
"""
function sanitize(y0,tsteps,tout)
    if ~isa(y0,AbstractVector)
        scalar = true
        y0 = [y0]
    else
        scalar = false
    end
    y0,tsteps,tout,scalar # TODO
end

"""
Converts the output to better output.

TODO: implement
"""
function desanitize(tout, yout, scalar=false)
    if scalar
        yout = squeeze(yout,1)
    end
    tout, yout
end

#################
# Julian mid level
################

"""
Solves the a problem `odep` and returns the solution (tout, yout),
"""
function solve(odep::ProblemFixedStep)
    @unpack odep.ode: y0, tspan
    T = typeof(tspan[1])
    dof = length(y0)
    N = length(odep)

    tout = zeros(T,N)
    tout[1] = tspan[1]
    yout = similar(y0,dof,N)
    yout[:,1] = y0
    for (status,st) in odep
        status::Status
        # write output
        tout[st.step] = st.t
        # allocates: yout[:,st.step] = st.y
        for i=1:dof
            yout[i,st.step] = st.y[i]
        end
    end
    return tout, yout
end
function solve(odep_d::Dense)
    @unpack odep_d.odep.ode: y0, tspan
    T = typeof(tspan[1])
    dof = length(y0)
    N = length(odep_d)
    # tout and yout are part of Dense
    for _ in odep_d
    end
    return odep_d.tout, odep_d.yout

end
function solve(odep::ProblemAdaptiveStep)
    @unpack odep.ode: y0, tspan
    tout = [tspan[1]]
    yout = copy(y0)
    for (status,st) in odep
        status::Status
        # write output
        push!(tout, st.t)
        append!(yout, st.y)
    end
    return tout, reshape(yout, length(y0), length(tout))
end

#################
# Julian high level
################

"""
Calls a solver to solve an ODE specified by an out-of-place F.
"""
ode{S<:ASolver}(Solver::Type{S}, F, y0, tsteps; opts...) = ode!(Solver, inplace(F), y0, tsteps; opts...)

# not sure about the naming...
function ode!{S<:ASolver}(Solver::Type{S}, F!, y0, tsteps; opts...)
    opts = Dict{Symbol,Any}(opts)
    get!(opts, :tsteps, tsteps)
    tout = get(opts, :tout, tsteps)
    y0,tsteps,tout,scalar = sanitize(y0,tsteps,tout)
    opts[:tsteps] = tsteps

    odep = Problem(ExplicitODE((tsteps[1],tsteps[end]), y0, F!),
                   Solver,
                   opts)
    odep_d = Dense(odep, tout)
    tout, yout = solve(odep_d)
    desanitize(tout, yout, scalar)
end


##################
# MATLAB-like interface
##################

ode4(F,y0,tsteps; opts...) = ode(RK4, F, y0, tsteps; opts...)
ode23s(F,y0,tsteps; opts...) = ode(Rosenbrock23s{:v2}, F, y0, tsteps; opts...)
