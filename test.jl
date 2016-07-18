using Base.Test

if !isdefined(:IVP)
    include("ODE.jl")
end
# Each entry is a test-set
# [F!, y0, tsteps, analytic_solution]
testsets = Vector[
                  Any[(t,y,dy)->dy[1] = 6.0, [0.], [0:.1:1;], t->6t],
                  Any[(t,y,dy)->dy[1]=2t, [0.], [0:.001:1;], t->t.^2],
                  Any[(t,y,dy)->dy[1]=y[1], [1.], [0:.001:1;], t->e.^t],
                  Any[(t,y,dy)->dy[1]=y[1], [1.], [1:-.001:0;], t->e.^(t-1)],
                  Any[(t,y,dy)-> (dy[1]= -y[2]; dy[2]= y[1]), [1., 2.], [0:.00001:2*pi;],
                         t->[cos(t)-2*sin(t) 2*cos(t)+sin(t)]' ]
                  ]
tol = 1e-3

# pick which solvers to use
fixedstep_solvers = Dict{Symbol,Any}()
fixedstep_solvers[:RK4] = RK4
fixedstep_solvers[:AdamMoulton] = AdamMoulton
fixedstep_solvers[:AdamBashforth] = AdamBashforth
merge!(fixedstep_solvers, rk_fixedstep_solvers)

@show solver_fixedstep = fixedstep_solvers[:RK4]
if solver_fixedstep in [AdamMoulton,AdamBashforth,RKFixedStep{:feuler}]
    tol = 1e-2
end
@show solver_adaptive = [Rosenbrock23s{:v1}, Rosenbrock23s{:v2}][2]

i=5
println("Fixed step stepper")
for i = 1:length(testsets)
    gc()
    println(i)
    F!,y0,tsteps,ana = testsets[i]
    tspan = (tsteps[1],tsteps[end])
    diffeq = ExplicitODE(tspan, y0, F!)
    opts = Dict(:tsteps=>tsteps, :order => 5)
    odep = Problem(diffeq,solver_fixedstep,opts)

    @inferred solve(odep)
# @show onestep_counter_RK4
# @show onestep_counter_rk
    @time tout, yout = solve(odep)

    if size(yout,1)==1
        yout = squeeze(yout,1)
    end
    @assert length(tout)==length(tsteps)
    @assert maximum(abs(yout-ana(tout))) < tol
    println("Error: $(maximum(abs(yout-ana(tout))))")
end

#error("asdf")

println("\nFixed step stepper with dense output")
for i = 1:length(testsets)
    gc()
    println(i)
    F!,y0,tsteps,ana = testsets[i]
    tspan = (tsteps[1],tsteps[end])
    diffeq = ExplicitODE(tspan, y0, F!)
    solver = solver_fixedstep
    opts = Dict(:tsteps=>tsteps,:order => 5)
    odep = Problem(diffeq,solver,opts)

    # output times
    tout = linspace(tspan[1],tspan[2], 17)
    dense = Dense(odep,tout)

    @inferred solve(dense)
    @time tout,yout = solve(dense)

    if size(yout,1)==1
        yout = squeeze(yout,1)
    end
    @assert length(tout)==17
    @assert maximum(abs(yout-ana(tout))) < tol
end

#error("asdf")


println("\nAdaptive step stepper")
i=3
for i = 1:length(testsets)
    solver = solver_adaptive
    gc()
    println(i)
    F!,y0,tsteps,ana = testsets[i]

    N = length(tsteps)
    dof = length(y0)
    tspan = (tsteps[1],tsteps[end])
    T = eltype(tspan)

    diffeq = ExplicitODE(tspan, y0, F!)
    opts = get_solver_defopts(diffeq,solver)
    opts[:reltol] = 1e-6

    odep = Problem(diffeq,solver,opts)

    @inferred solve(odep)
    @time tout, yout = solve(odep)

    if size(yout,1)==1
        yout = squeeze(yout,1)
    end

    @assert maximum(abs(yout-ana(tout))) < tol
end


println("\nAdaptive step stepper with dense output")
for i = 1:length(testsets)
    solver = solver_adaptive
    gc()
    println(i)
    F!,y0,tsteps,ana = testsets[i]
    tspan = (tsteps[1],tsteps[end])

    diffeq = ExplicitODE(tspan, y0, F!)
    opts = get_solver_defopts(diffeq,solver)
    opts[:reltol] = 1e-6

    odep = Problem(diffeq,solver,opts)

    # output times
    tout = linspace(tspan[1],tspan[2], 17)
    dense = Dense(odep, tout)

    @inferred solve(dense)
    @time tout,yout = solve(dense)

    if size(yout,1)==1
        yout = squeeze(yout,1)
    end
    @assert length(tout)==17
    @assert maximum(abs(yout-ana(tout))) < tol
end


use_scalar_fns = false

println("\nHigh level interface. Using with out-of-place functions")

if use_scalar_fns
    # Note that ode4_old is faster than the new implementation when
    # using straight floats.  Even faster than the in-place functions
    # above.  This is the price we pay for making them all vectors
    # internally.  However, I suspect it's a small price to pay in all
    # but crazy examples.
    testsets = Vector[
                      Any[(t,y)->6.0, 0., [0:.1:1;], t->6t],
                      Any[(t,y)->2t, 0., [0:.001:1;], t->t.^2],
                      Any[(t,y)->y, 1., [0:.001:1;], t->e.^t],
                      Any[(t,y)->y, 1., [1:-.001:0;], t->e.^(t-1)],
                      Any[(t,y)->[-y[2]; y[1]], [1., 2.], [0:.00001:2*pi;], # this is always non-scalar
                          t->[cos(t)-2*sin(t) 2*cos(t)+sin(t)]' ]
                      ]
else
    testsets = Vector[
                      Any[(t,y)->[6.0], [0.], [0:.1:1;], t->(6t)'],
                      Any[(t,y)->[2t], [0.], [0:.001:1;], t->(t.^2)'],
                      Any[(t,y)->y, [1.], [0:.001:1;], t->(e.^t)'],
                      Any[(t,y)->y, [1.], [1:-.001:0;], t->(e.^(t-1))'],
                      Any[(t,y)->[-y[2]; y[1]], [1., 2.], [0:.00001:2*pi;],
                          t->[cos(t)-2*sin(t) 2*cos(t)+sin(t)]' ]
                      ]
end

#ODE4  Solve non-stiff differential equations, fourth order
#   fixed-step Runge-Kutta method.
#
#   [T,X] = ODE4(ODEFUN, X0, TSPAN) with TSPAN = [T0:H:TFINAL]
#   integrates the system of differential equations x' = f(t,x) from time
#   T0 to TFINAL in steps of H with initial conditions X0. Function
#   ODEFUN(T,X) must return a column vector corresponding to f(t,x). Each
#   row in the solution array X corresponds to a time returned in the
#   column vector T.
function ode4_old(F, x0, tspan)
    h = diff(tspan)
    x = Array(typeof(x0), length(tspan))
    x[1] = x0

    midxdot = Array(typeof(x0), 4)
    for i = 1:length(tspan)-1
        # Compute midstep derivatives
        midxdot[1] = F(tspan[i],         x[i])
        midxdot[2] = 2*F(tspan[i]+h[i]./2, x[i] + midxdot[1].*h[i]./2)
        midxdot[3] = 2*F(tspan[i]+h[i]./2, x[i] + midxdot[2].*h[i]./2)
        midxdot[4] = F(tspan[i]+h[i],    x[i] + midxdot[3].*h[i])

        # Integrate
        x[i+1] = x[i] + 1/6 .*h[i].*sum(midxdot)
    end
    return tspan, x
end

import ODE

for i = 1:length(testsets)
    gc()
    println(i)
    F,y0,tsteps,ana = testsets[i]

    # output times
    tout = linspace(tsteps[1],tsteps[end], 17)

    ode(RK4, F,y0,tsteps)
    gc()
    print("RK4:        ")
    @time tout,yout = ode(RK4, F, y0, tsteps)
    @assert length(tout)==length(tsteps)
    @assert maximum(abs(yout-ana(tout))) < tol

    ode4(F,y0,tsteps)
    gc()
    print("ode4:       ")
    @time tout,yout = ode4(F,y0,tsteps)
    @assert length(tout)==length(tsteps)
    @assert maximum(abs(yout-ana(tout))) < tol

    ode4_old(F,y0,tsteps)
    gc()
    print("ODE.ode4:   ")
    @time tout,yout = ode4_old(F,y0,tsteps)
    @assert length(tout)==length(tsteps)
    if eltype(yout)<:Array
        yout = hcat(yout...)
    end
    @assert maximum(abs(yout-ana(tout))) < tol

    ode23s(F,y0,tsteps, reltol=1e-6)
    gc()
    print("ode23s:     ")
    @time tout,yout = ode23s(F,y0,tsteps, reltol=1e-6)
    @assert length(tout)==length(tsteps)
    @assert maximum(abs(yout-ana(tout))) < tol

    # Test against ode23s from ODE.jl
    ODE.ode23s(F,y0,tsteps)
    gc()
    print("ODE.ode23s: ")
    @time tout,yout = ODE.ode23s(F,y0,tsteps, reltol=1e-6)
    yout = hcat(yout...)
    @assert maximum(abs(yout-ana(tout))) < tol
end
