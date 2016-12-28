abstract ODEIterAlgorithm <: AbstractODEAlgorithm
immutable feuler <: ODEIterAlgorithm end
immutable rk23 <: ODEIterAlgorithm end
immutable rk45 <: ODEIterAlgorithm end
immutable feh78 <: ODEIterAlgorithm end
immutable ModifiedRosenbrock <: ODEIterAlgorithm end
immutable midpoint <: ODEIterAlgorithm end
immutable heun <: ODEIterAlgorithm end
immutable rk4 <: ODEIterAlgorithm end
immutable feh45 <: ODEIterAlgorithm end

function solve{uType,tType,isinplace,algType<:ODEIterAlgorithm,F}(
      prob::AbstractODEProblem{uType,tType,isinplace,F},
      alg::algType,timeseries=[],ts=[],ks=[];dense=true,
      timeseries_errors=true,dense_errors=false,kwargs...)

    tspan = prob.tspan

    o = Dict{Symbol,Any}(kwargs)
    u0 = prob.u0
    o[:T] = tspan[end]

    if typeof(u0) <: Number
        u = [u0]
    else
        u = deepcopy(u0)
    end

    sizeu = size(u)

    opts = buildOptions(o,ODEJL_OPTION_LIST,ODEJL_ALIASES,ODEJL_ALIASES_REVERSED)

    if !isinplace && typeof(u)<:AbstractArray
        f! = (t,u,du) -> (du[:] = prob.f(t,u))
    else
        f! = prob.f
    end
    ode  = ODE.ExplicitODE(tspan[1],u,f!)
    # adaptive==true ? FoA=:adaptive : FoA=:fixed #Currently limied to only adaptive
    FoA = :adaptive
    if algType <: rk23
        solver = ODE.RKIntegrator{FoA,:rk23}
    elseif algType <: rk45
        solver = ODE.RKIntegrator{FoA,:dopri5}
    elseif algType <: feh78
        solver = ODE.RKIntegrator{FoA,:feh78}
    elseif algType <: ModifiedRosenbrock
        solver = ODE.ModifiedRosenbrockIntegrator
    elseif algType <: feuler
        solver = ODE.RKIntegratorFixed{:feuler}
    elseif algType <: midpoint
        solver = ODE.RKIntegratorFixed{:midpoint}
    elseif algType <: heun
        solver = ODE.RKIntegratorFixed{:heun}
    elseif algType <: rk4
        solver = ODE.RKIntegratorFixed{:rk4}
    elseif algType <: feh45
        solver = ODE.RKIntegrator{FoA,:rk45}
    end
    out = ODE.solve(ode;solver=solver,opts...)
    y = out.y
    t = out.t
    dy = out.dy
    if length(out.y[1])==1
        tmp = Vector{eltype(out.y[1])}(length(out.y))
        tmp_dy = Vector{eltype(out.dy[1])}(length(out.dy))
        for i in 1:length(out.y)
            tmp[i] = out.y[i][1]
            tmp_dy[i] = out.dy[i][1]
        end
        y = tmp
        dy = tmp_dy
    end

    #saveat_idxs = find((x)-> x ∈ saveat,ts)
    #t_nosaveat = view(ts,symdiff(1:length(ts),saveat_idxs))
    #u_nosaveat = view(timeseries,symdiff(1:length(ts),saveat_idxs))

    if dense
        interp = (tvals) -> common_interpolation(tvals,t,y,dy,alg,f!)
    else
        interp = (tvals) -> nothing
    end

    build_solution(prob,alg,t,y,
                      dense=dense,k=dy,interp=interp,
                      timeseries_errors = timeseries_errors,
                      dense_errors = dense_errors)
end

const ODEJL_OPTION_LIST = Set([:tout,:tstop,:reltol,:abstol,:minstep,:maxstep,:initstep,:norm,:maxiters,:isoutofdomain])
const ODEJL_ALIASES = Dict{Symbol,Symbol}(:minstep=>:dtmin,:maxstep=>:dtmax,:initstep=>:dt,:tstop=>:T,:maxiters=>:maxiters)
const ODEJL_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEJL_ALIASES])

function buildOptions(o,optionlist,aliases,aliases_reversed)
    dict1 = Dict{Symbol,Any}([Pair(k,o[k]) for k in (keys(o) ∩ optionlist)])
    dict2 = Dict([Pair(aliases_reversed[k],o[k]) for k in (keys(o) ∩ values(aliases))])
    merge(dict1,dict2)
end

"""
common_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function common_interpolation(tvals,ts,timeseries,ks,alg,f)
  idx = sortperm(tvals)
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  vals = Vector{eltype(timeseries)}(length(tvals))
  for j in idx
    t = tvals[j]
    i = findfirst((x)->x>=t,ts[i:end])+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[i] == t
      vals[j] = timeseries[i]
    elseif ts[i-1] == t # Can happen if it's the first value!
      vals[j] = timeseries[i-1]
    else
      dt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/dt
      vals[j] = common_interpolant(Θ,dt,timeseries[i-1],timeseries[i],ks[i-1],ks[i],alg)
    end
  end
  vals
end

"""
common_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function common_interpolation(tval::Number,ts,timeseries,ks,alg,f)
  i = findfirst((x)->x>=tval,ts) # It's in the interval ts[i-1] to ts[i]
  if ts[i] == tval
    val = timeseries[i]
  elseif ts[i-1] == tval # Can happen if it's the first value!
    push!(vals,timeseries[i-1])
  else
    dt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/dt
    val = common_interpolant(Θ,dt,timeseries[i-1],timeseries[i],ks[i-1],ks[i],alg)
  end
  val
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190
"""
function common_interpolant(Θ,dt,y₀,y₁,k₀,k₁,alg) # Default interpolant is Hermite
  (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k₀ + Θ*dt*k₁)
end

export ODEIterAlgorithm, feuler, rk23, feh45, feh78, ModifiedRosenbrock,
      midpoint, heun, rk4, rk45
