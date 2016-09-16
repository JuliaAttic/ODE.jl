var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "Pages = [\n    \"tutorials/euler_integrator.md\",\n    \"man/basics.md\",\n    \"man/base.md\"\n    ]"
},

{
    "location": "index.html#ODE.jl-1",
    "page": "Home",
    "title": "ODE.jl",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Top-level-interface-1",
    "page": "Home",
    "title": "Top level interface",
    "category": "section",
    "text": "If you are looking to getting a solution without the additional hussle of handling an iterator we provide the wrappers ODE.odeXX.  They provide a simplistic way of handling explicit differential equations of the form y'=F(t,y) with y being either a Number or an AbstractArray of any dimension.  Below we solve a simple initial value problem given by y'=y with initial data at t0=0.0 and y0=1.0 on the interval [t0,1].using ODE\ntspan  = [0.0,1.0]\ny0     = 1.0\nF(t,y) = y\n(t,y)  = ODE.ode(F,y0,tspan)The vectors t and y store the time and solution values at the corresponding times."
},

{
    "location": "index.html#Solve-interface-1",
    "page": "Home",
    "title": "Solve interface",
    "category": "section",
    "text": "ODE.ode only supports explicit differential equations defined as y'=F(t,y), for more advenced uses consider using ODE.solve, which was designed to work with a variety of other types of initial value problems and is optimized for better performance.  First we have to define an initial value problem, in our case this is an explicit differential equation y'=y with inital data y0=[1.0] given at the time t0=0.0.using ODE\nt0  = 0.0\ny0  = [1.0]\nF!(t,y,dy) = dy[1]=y[1]\node = ODE.ExplicitODE(t0,y0,F!)Note that unlike in ODE.ode we now have to supply an in place function F! instead of an explicit function F.  We can solve the ODE problem ode by simply callingsol = ODE.solve(ode, tstop = 1)This returns a Solution type, which stores the solution.  You probably noticed that we passed a keyword argument tstop, this is the final time of integration which we have to specify because tstop defaults to Inf and the integration would carry on forever.  You can access the solution with(t,y) = sol.t, sol.yYou can change the default algorithm (Runge-Kutta (4,5)) by passing an optional argument solversol = ODE.solve(ode, tstop = 1, solver = ODE.RKIntegratorAdaptive{:dopri5})For other options accepted by solve see Options below.You might still find this interface limiting.  First of all, it stores all the results, so if you are only interested in the final value of y it still stores all the intermediate steps.  Secondly, you cannot process the results on the fly (e.g. plot the current state of a solution).  If you need more control you should consider using the iterator interface."
},

{
    "location": "index.html#Iterator-interface-1",
    "page": "Home",
    "title": "Iterator interface",
    "category": "section",
    "text": "To offeset the limitations of the ODE.ode interface we implemented a general.  We use the same problem as before as an exampleusing ODE\nt0  = 0.0\ny0  = [1.0]\nF!(t,y,dy) = dy[1]=y[1]\node = ODE.ExplicitODE(t0,y0,F!)Now we have full flow control over the solver, we can analyze the intermediate results or interrupt the integration at any point.for (t,y) in ODE.iterate(ode)\n    @show (t,y)\n    if t > 1\n        break\n    end\nendNote that we had to break the loop because sol would keep producing the results.  To set the final integration time and other parameters of the integrator integ we can pass optional arguments to ODE.solver.for (t,y) in ODE.iterate(ode; tstop = 1)\n    @show (t,y)\nendThis approach has the added benefit of the solution never exceeding the final time.  Both ODE.iterate and ODE.solve support the same options, so you can easily change the method of integration with the keyword solver.Apart from the time and value (t,y) the ODE.solve also returns the derivative, you can retrive it as the third argument in the returned tuple.  In the following example we use it compute the absolute residual error (zero in this case).for (t,y,dy) in ODE.iterate(ode; tstop = 1)\n    err = norm(y-dy)\n    @show err\nendWith tstop specified we can also get all results at once using collect and other constructs working on iterators (e.g. generators).  For exampleiter = ODE.iterate(ode; tstop = 1)\nsolution = collect(iter)returns a vector of triples (t,y,dy).  Or if you only wan the first component of a solution you could simply usey1 = collect(y[1] for (t,y) in iter)There are, however, several caveats that you should take into account:Each time the iterator is collected the differential equation is actually solved, which has potentially high computational cost and might be inefficient.\nThe length is undefined for the result of ODE.iterate, because we don't know a priori how many steps the integration will require (especially in the case of adaptive solvers).  This means that the functions requireing length might not work.  For the same reason there are no getindex methods."
},

{
    "location": "index.html#Options-1",
    "page": "Home",
    "title": "Options",
    "category": "section",
    "text": "Both ODE.ode and ODE.solve accept the following keyword arguments.integ: the type of integrator to use, defaults to a adaptive Runge-Kutta method of order 4/5.  To see the list of available integrators see Integrators.\ninitstep: The initial step size, defaults to eps(T)^(1/3).\ntstop: The final integration time, never exceeded by the integrator.  In case of ODE.ode(F,y0,tspan) this option defaults to the last element of tspan if it is a vector.  In ODE.solve the default is tstop=Inf.  If tstop is smaller then t0 the integration runs backwards in time.Apart from these general options, each integrator has its own keyword arguments.  In particular all integrators with adaptive step size can be cotrolled withreltol, abstol: The relative and absolute error tolerances.  The solution guarantees that at each step we have norm((y-yc)*reltol.+abstol)<=1, where yc is a true solution to and IVP.  Defaults are reltol=eps(T)^(1/3)/10, abstol=eps(T)^(1/2)/10.\nnorm: The norm used to measure error in the formula above, defaults to y->Base.vecnorm(y,Inf).  You can specify it to assign different weights to different components of y.\nminstep, maxstep: Minimal and maximal stepsize for the integrator.  If at any point the stepsize exceeds these limits the integrator will yield an error and cease producing results.  Deafaults are minstep=10*eps(T) and maxstep=1/minstep.\nmaxiters: The number of iterations before the integrator ceases to work, defaults to Inf.  Useful as a safeguard from iterator continuing ad infinitum.\nisoutofdomain: Applied to each component of y, if isoutofdomain(y[i])==true the integrator stops.  Defaults to Base.isnan.Apart from these, each integrator may support additional options."
},

{
    "location": "index.html#Integrators-1",
    "page": "Home",
    "title": "Integrators",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Explicit-Runge-Kutta-integrators-1",
    "page": "Home",
    "title": "Explicit Runge-Kutta integrators",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Rosenbrock-methods-1",
    "page": "Home",
    "title": "Rosenbrock methods",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Backwards-differential-formula-(BDF)-methods-1",
    "page": "Home",
    "title": "Backwards differential formula (BDF) methods",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#???-1",
    "page": "Home",
    "title": "???",
    "category": "section",
    "text": ""
},

{
    "location": "man/basics.html#",
    "page": "Basics",
    "title": "Basics",
    "category": "page",
    "text": ""
},

{
    "location": "man/basics.html#Basic-usage-1",
    "page": "Basics",
    "title": "Basic usage",
    "category": "section",
    "text": "Consider an ODE y=y"
},

{
    "location": "man/base.html#",
    "page": "Base",
    "title": "Base",
    "category": "page",
    "text": "CurrentModule = ODE"
},

{
    "location": "man/base.html#Base-1",
    "page": "Base",
    "title": "Base",
    "category": "section",
    "text": "The file base.jl implements the most basic iterator infrastructure for solvers and the definitions of the types representing general IVP (initial value problem) and solvers."
},

{
    "location": "man/base.html#ODE.solve",
    "page": "Base",
    "title": "ODE.solve",
    "category": "Function",
    "text": "solve(ivp::IVP; solver=RKIntegratorAdaptive{:rk45}, opts...)\n\nSolve the initial value problem ivp using an algorithm solver (defaults to Runge-Kutta (4,5) integrator).  One can pass additional options to the solver via keyword arguments to solve (here denoted as options).  The output is a Solution type (currently simply a tuple of vectors (Vector{T},Vector{Y}), where T,Y=eltype(ivp)).\n\n\n\n"
},

{
    "location": "man/base.html#ODE.iterate",
    "page": "Base",
    "title": "ODE.iterate",
    "category": "Function",
    "text": "iterate(ivp::IVP; solver=RKIntegratorAdaptive{:rk45}, opts...)\n\nIterate creates an iterable Problem instance from an IVP instance (specifying the math) and from a Type{AbstractSolver} (the numerical integrator).  The simplest use case is\n\nfor (t,y,dy) in iterate(...)\n    # do something with t, y an dy\nend\n\nIf the integration interval, defined by the keyword argument tstop, is finite you can request all the results at once by calling\n\ncollect(iterate(...)) # => Vector{Tuple{T,Y,Y}}\n\nNotes:\n\nusually solvers require the ivp to be in a certain form, say an\n\nExplicitODE. - the second argument is the Type of the solver and not an instance.   The instance of the solve can only be created together with the   ivp as their type parameters need to match.\n\nInput:\n\nivp::IVP\nS::Type{AbstractSolver}\n\nOutput:\n\n::Problem\n\n\n\n"
},

{
    "location": "man/base.html#General-functions-for-solving-initial-value-problems-1",
    "page": "Base",
    "title": "General functions for solving initial value problems",
    "category": "section",
    "text": "solve\niterate"
},

{
    "location": "man/base.html#ODE.AbstractIVP",
    "page": "Base",
    "title": "ODE.AbstractIVP",
    "category": "Type",
    "text": "AbstractIVP{T,Y}\n\nThe abstract supertype of all IVPs (initial value problems).  The type parameters T and Y correspond to the types of time and state variable respectively.\n\n\n\n"
},

{
    "location": "man/base.html#ODE.IVP",
    "page": "Base",
    "title": "ODE.IVP",
    "category": "Type",
    "text": "IVP{T,Y,F,G,J} <: AbstractIVP{T,Y}\n\nDefines the mathematical part of an IVP (initial value problem) specified in the general form:\n\nF(t, y) =  G(t, y, dy) with y(t0)= y0\n\nDepending on the combination of the parameters this type can represent a wide range of problems, including ODE, DAE and IMEX.  Nevertheless not all solvers will support any combinations of F and G.  Note that not specifying G amounts to G=dy/dt.\n\ntspan – tuple (start_t,end_t)\ny0 – initial condition\nF! – in-place F function F!(t,y,res).  If F=0 set to nothing.\nG! – in-place G function G!(t,y,dy,res).  If G=dy/dt then         set to nothing (or dy if the solver supports this).  Can         also be a mass matrix for a RHS M dy/dt\nJ! – in-place Jacobian function J!(t,y,dy,res).\n\nTODO: how to fit the sparsity pattern in J?\n\n\n\n"
},

{
    "location": "man/base.html#ODE.ExplicitODE",
    "page": "Base",
    "title": "ODE.ExplicitODE",
    "category": "Constant",
    "text": "typealias ExplicitODE{T,Y,F,J} IVP{T,Y,F,Void,J}\n\nCan be constructed by calling\n\nODE.ExplicitODE(t0,y0,F!;J!=jacobian))\n\nExplicit ODE representing the problem\n\ndy = F(t,y) with y(t0)=y0\n\nt0, y0: initial conditions\nF!: in place version of F called by F!(t,y,dy)\nJ!: (optional) computes J=dF/dy in place, called with J!(t,y,J)\n\n\n\n"
},

{
    "location": "man/base.html#ODE.ImplicitODE",
    "page": "Base",
    "title": "ODE.ImplicitODE",
    "category": "Constant",
    "text": "Implicit ODE representing the problem\n\nG(t,y,dy)=0 with y(t0)=y0 and optionally y'(t0)=dy0\n\nt0, y0: initial conditions\nG!: in place version of G called by G!(res,t,y,dy),     returns residual in-place in res.\nJ!: (optional) computes J=dF/dy+a*dF/dy' for prescribed a, called with J!(out,t,y,dy,a).     Returns Jacobian in-place in out.\n\n\n\n"
},

{
    "location": "man/base.html#Predefined-types-of-initial-value-problems-1",
    "page": "Base",
    "title": "Predefined types of initial value problems",
    "category": "section",
    "text": "AbstractIVP\nIVP\nExplicitODE\nImplicitODE"
},

{
    "location": "man/base.html#ODE.AbstractSolver",
    "page": "Base",
    "title": "ODE.AbstractSolver",
    "category": "Type",
    "text": "AbstractSolver\n\nThe supertype of anything which can get you to a solution of a IVP. Subtypes include: AbstractIntegrators but also DenseOutput\n\n\n\n"
},

{
    "location": "man/base.html#ODE.AbstractIntegrator",
    "page": "Base",
    "title": "ODE.AbstractIntegrator",
    "category": "Type",
    "text": "The abstract type of the actual algorithm to solve an IVP.\n\n\n\n"
},

{
    "location": "man/base.html#ODE.AbstractState",
    "page": "Base",
    "title": "ODE.AbstractState",
    "category": "Type",
    "text": "AbstractState keeps the temporary data (state) for the iterator Problem{::AbstractIntegrator}.\n\n\n\n"
},

{
    "location": "man/base.html#ODE.Solution",
    "page": "Base",
    "title": "ODE.Solution",
    "category": "Type",
    "text": "Stores a solution to the ivp.\n\n\n\n"
},

{
    "location": "man/base.html#Solver-architecture-1",
    "page": "Base",
    "title": "Solver architecture",
    "category": "section",
    "text": "AbstractSolver\nAbstractIntegrator\nAbstractState\nSolutionThe fallback constructor for AbstractSolver(ivp::IVP;opts...) ensures that an error is thrown if a solver is constructed for an unsupported type of the given IVP."
},

{
    "location": "man/base.html#Fallback-functions-for-solvers-1",
    "page": "Base",
    "title": "Fallback functions for solvers",
    "category": "section",
    "text": "Base.length(::AbstractSolver)\noutput(::AbstractState)\nBase.eltype{T,Y}(::Type{AbstractIVP{T,Y}})"
},

]}
