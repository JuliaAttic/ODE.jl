# Shows the low-level iteration

plot_yes = false

if !isdefined(:IVP)
    include("ODE.jl")
end

F!,y0,tsteps,ana = Any[(t,y,dy)-> (dy[1]= -y[2]; dy[2]= y[1]), [1., 2.], [0:.00001:2*pi;],
                         t->[cos(t)-2*sin(t) 2*cos(t)+sin(t)]' ]
tol = 1e-3

tspan = (tsteps[1],tsteps[end]) # tuple (t0, t_end)

# construct a IVP instance from above functions:
diffeq = ExplicitODE(tspan, y0, F!)

# pick a solver
solver = [RK4, Rosenbrock23s{:v1}, Rosenbrock23s{:v2}][3]

# Get default options for a IVP + Solver combination
opts = get_solver_defopts(diffeq,solver)
# Unknown options will be ignored by a solver.
opts[:reltol] = 1e-6
opts[:tsteps] = tsteps

# Create the `Problem` from the equation, solver + options.
#
# `opts` can be left away.
odep = Problem(diffeq,solver,opts)


### Dense with events output:

# output times
tout_ = linspace(tspan[1],tspan[2], 13)

ef1(t,y) = y[1]-0.5
ef2(t,y) = y[1] + 0.2*t -3
isterminal = [false, true]

# Just wrap the `Problem` with Dense
events = DenseEvents(odep,tout_, (ef1,ef2), isterminal=isterminal)
dof = length(y0)
N = length(events)


etout = Float64[]
eyout = similar(y0,0)
for (i,eid,t,y) in events
    if eid!=0
        push!(etout, t)
        append!(eyout, y)
    end
end
eyout = reshape(eyout, dof, length(etout))

yout = events.dense.yout
tout = events.dense.tout
@assert maximum(abs(eyout-ana(etout))) < tol
@assert maximum(abs(yout-ana(tout))) < tol

if plot_yes && VERSION<v"0.5-"
    eval(:(using Plots))
    mp = plot(tout, vec(yout[1,:]), marker=:hex, label="sol")
    plot!(mp, tout, -Float64[ef1(tout[i],yout[:,i])-yout[1,i] for i=1:length(tout)], label="ef1")
    plot!(mp, tout, -Float64[ef2(tout[i],yout[:,i])-yout[1,i] for i=1:length(tout)], label="ef2")
    scatter!(mp, etout, vec(eyout[1,:]), color=:green, label="events")
end
