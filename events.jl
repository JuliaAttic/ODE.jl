# An event happen when a function f(t,y)==0.  Thus this is
# equivalent of root finding.

# TODO:
# - this duplicates much of dense.jl: DenseEvents is a superset of
#   Dense.  Can this be refactored? --> if there is the :all option, then yes.
# - allow functions of type f(t,y,dy)?  -> that makes interpolation harder.


"Direction of zero crossings for event functions (roots)."
@enum Dir neg2pos=1 pos2neg=-1 both=0

"""
To get dense output & events, wrap a `Problem` with `DenseEvents`.

# Event functions
Have a signature `f(t,y,dy) --> Float64`, an event occurs at roots of
f.  The type takes a tuple of event functions.  Further options:

- `isterminal::BitVector` -- if true that event terminates the
                             integration (default `false`)
- `direction::Dir` -- whether only one or both directions of zero
                      crossing trigger an event (default `both`)

# Constructor
```
DenseEvents{OP<:Problem}(odep::OP, tout::AbstractVector, efs=();
                            efs_abstol=ones(length(efs))*1e-6,
                            isterminal=falses(length(efs)),
                            direction=[both for i=1:length(efs)])
```

TODO: remove tout, yout from this type?
"""
type DenseEvents{OP<:Problem, T, Y2<:AbstractMatrix, EF} # TODO: make immutable?
    dense::Dense{OP,T,Y2}
    points::Symbol          # to store original value of ds.points, as it will be overwritten
    # events
    efs::EF         # tuple containing the eventfunctions
    efs_abstol::Vector{Float64} # tolerance with which to determine root
    isterminal::Vector{Bool} # terminate integration if set
    direction::Vector{Dir}   # which dir of zero-crossings
end
# TODO:
#  - get tout from opts?
#  - think about wrapping solver which only have dense output.  Is that possible?
#  - add :all option
function DenseEvents{OP<:Problem,T}(odep::OP, tout::AbstractVector{T},
                                  efs=(), opts::Dict{Symbol}=Dict{Symbol,Any}();
                                  efs_abstol=ones(length(efs))*1e-6,
                                  isterminal=falses(length(efs)),
                                  direction=[both for i=1:length(efs)])
    dense = Dense(odep, tout, opts)
    points = dense.points
    dense.points = :all # needed in nested iteration
    if any(direction.!=both)
        error("not implemented") # TODO
    end
    DenseEvents{OP,T,typeof(dense.yout),typeof(efs)}(dense,points,efs,efs_abstol,isterminal,direction)
end

Base.length(events::DenseEvents) = length(events.dense.tout)

type DenseEventsState{ST<:AState,T,Y2}
    ds::DenseState{ST}
    efvals::Vector{Float64}  # values of event functions
    efid::Vector{Int}       # indices of triggered events
    eft::Vector{T}          # times of triggered events
    terminate::Bool        # terminate iteration
    ytmp::Y2                # to hold y at event.  2D array for type-stability
    t_tmp::T                # time at last step
end

function DenseEventsState{OP,T,Y2}(events::DenseEvents{OP,T,Y2})
    @unpack events: efs
    @unpack odep.ode: y0
    ds = DenseState(events.dense)
    st = ds.st
    ST = typeof(st)

    efvals = Array(Float64,length(efs))
    t0 = gett(st)
    y0 = gety(st)
    dy0 = getdy(st)
    for (i,ef) in enumerate(efs)
        efvals[i] = ef(t0,y0)
        # TODO: what if there is a root at t0?
    end
    efid = Int[]; sizehint!(efid, length(efs))
    eft = T[]; sizehint!(eft, length(efs))
    ytmp = similar(y0,length(y0),1)
    terminate = false
    t_tmp = t0

    DenseEventsState{ST,T,Y2}(ds, efvals, efid, eft, terminate, ytmp, t_tmp)
end


"""
Iteration item is:

- index into dense.tout, zero if a points==:all point
- index into event function tuple, zero if no event occured
- tout
- yout
"""
Base.start(events::DenseEvents) = DenseEventsState(events)
function Base.next(events::DenseEvents, es::DenseEventsState)
    @unpack events: dense, efs, efs_abstol, isterminal
    @unpack dense: yout, tout
    @unpack es: ds, efvals, eft, efid, ytmp, t_tmp
    @unpack ds: i, st
    @unpack events.dense.odep.ode: tspan

    tdir = sign(tspan[2]-tspan[1])

    # Iterate over `Dense`:
    #  - return dense values, unless
    #  - an event function changes sign

    # If no events to process, do a step with the dense stepper:
    while length(efid)==0
        # do a step with dense stepper (note that points=:all)
        (i,t,y), ds = next(dense, ds) # mutates ds and st

        # go through event functions to check for sign-change
        for (ii,ef) in enumerate(efs)
            v = ef(t,y)
            if sign(v)!=sign(efvals[ii])
                # TODO direction
                push!(efid,ii)
                # find zero
                rng = (t_tmp, t)
                push!(eft, findroot(rootfn, rng, efs_abstol[ii], ef, ytmp, st))
            end
            efvals[ii] = v
        end
        # return dense state if no events were found:
        if length(efid)==0
            if i==0 # points==:all point
                if events.points==:all
                    return (0, 0, t, y), es
                else
                    # do one more step with dense iterator
                end
            else
                return (i, 0, t, y), es
            end
        else # if events found then sort them, first is last
            perm = sortperm(eft,rev=true)
            eft[:] = eft[perm]
            efid[:] = efid[perm]
        end
        es.t_tmp = t
    end
    # if we get to here, then return events
    tt = pop!(eft)
    id = pop!(efid)
    es.terminate = isterminal[id]
    interpolate!(st, tt, ytmp)

    return  (0, id, tt, view(ytmp, :, 1)), es # return view for type-stability
    # TODO: return Status
end
function Base.done(events::DenseEvents, es::DenseEventsState)
    if done(events.dense, es.ds)
        return true
    elseif es.terminate # early termination
        # need to tidy up dense.yout and dense.tout as not all were populated
        #
        # TODO: this is ugly, also means that DenseEvents cannot be immutable
        @unpack events.dense: yout, tout
        @unpack es.ds: i
        i -= 1
        dof = size(yout,1)
        yout = reshape(yout, length(yout))
        resize!(yout, i*dof)
        events.dense.yout = reshape(yout, dof, i)
        resize!(tout, i)
        return true
    else
        return false
    end
end



# Used in the root finding algo:
rootfn(t, ef, ytmp, st) = (interpolate!(st,t,ytmp); ef(t, ytmp))
# Rootfinding from @pwl's PR49
"""
A simple bisection algorithm for finding a root of a solution f(x)=0
starting within the range x∈rng, the result is a point x₀ which is
located within the distance abstol from the true root of f(x)=0.  For
this algorithm to work we need f(rng[1]) to have a different sign then
f(rng[2]).
"""
function findroot(f,rng,abstol, args...)
    xl, xr = rng
    fl, fr = f(xl, args...), f(xr, args...)

    if fl*fr > 0 || xl > xr
        error("Inconsistent bracket")
    end

    while xr-xl > abstol
        xm = (xl+xr)/2
        fm = f(xm, args...)

        if fm*fr > 0
            xr = xm
            fr = fm
        else
            xl = xm
            fl = fm
        end
    end

    return (xr+xl)/2
end
