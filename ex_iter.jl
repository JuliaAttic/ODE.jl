# Shows the low-level iteration
include("ODE.jl")

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

### Just output the steps which were taken:
# (note if performance is important this should be inside a function)

tout = [tspan[1]]
yout = similar(y0)
yout[:] = y0
for (status,st) in odep
    # status is empty for now

    # write output
    push!(tout, gett(st))
    append!(yout, gety(st))
end
yout = reshape(yout, length(y0), length(tout))

# check against analytic solution
@assert maximum(abs(yout-ana(tout))) < tol

### Dense output:

# output times
tout = linspace(tspan[1],tspan[2], 17)

# Just wrap the `Problem` with Dense
odep_d = Dense(odep,tout)
dof = length(y0)
N = length(odep_d)


for (i,t,y) in odep_d #
    # Actually, no need to do anything as yout is part of odep_d.
    # Although, I'm not sure that's the best.
end
yout = odep_d.yout
@assert maximum(abs(yout-ana(tout))) < tol
