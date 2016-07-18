# Make dense output

"""

Note, only "specified" points get stored in the internal datastores.
"""
type Dense{OP<:Problem, T, Y2<:AbstractMatrix}
    odep::OP
    tout::Vector{T}
    yout::Y2
    # options
    points::Symbol  # [:all, :specified][2]
end
# TODO:
#  - get tout from opts?
#  - think about wrapping solver which only have dense output.  Is that possible?
function Dense{OP<:Problem}(odep::OP, tout::AbstractVector, opts::Dict{Symbol}=Dict{Symbol,Any}())
    merge!(odep.opts, opts)
    y0 = odep.ode.y0
    N = length(tout)
    dof = length(y0)
    T = eltype(odep.ode.tspan)
    yout = similar(y0,dof,N)
    points = get!(opts, :points, :specified)
    @assert points in [:all, :specified]
    Dense{OP,T,typeof(yout)}(odep,tout,yout, points)
end
Base.length(d::Dense) = length(d.tout)

type DenseState{ST<:AState}
    st::ST
    i::Int # index into tout
end
function DenseState(dense::Dense)
    st = init(dense.odep)
    DenseState{typeof(st)}(st, 1)
end

# Iteration item is:
#  - index into dense.tout, zero if a points==:all point
#  - tout
#  - yout
Base.start(dense::Dense) = DenseState(dense)
function Base.next(dense::Dense, ds::DenseState)
    @unpack dense: yout, points
    tspan = dense.odep.ode.tspan
    tdir = sign(tspan[2]-tspan[1])
    i = ds.i
    tout = dense.tout[i]


    st = ds.st
    t = gett(st)
    while tout*tdir>=t*tdir && !done(dense.odep, st)
        next(dense.odep, st)
        t = gett(st)
        if points==:all
            return (0, t, gety(st)), ds
        end
    end
    # now we got a st such that tout is included
    interpolate!(st, tout, yout, i)


    ds.i += 1
    return (ds.i, tout, view(yout, :, i)), ds
    # TODO: return Status
end
Base.done(dense::Dense, ds::DenseState) = ds.i==length(dense.tout)+1

"""
Interpolates the solution encoded in the state to the time `tout` and
writes the result into `yout`.  The default implementation uses a 3rd
order Hermite interpolation.  If desired, this function can be
overloaded for a particular solver by dispatching on its state.

Input:

- `st::AState` -- The state of the ode iteration.  It has to hold all
                  the information needed for the interpolation.  Thus
                  a state type has to be designed with this in mind.
- `tout` -- time of requested output.  Note that
            `gett_pre(st)<=tout<=gett(st)`!
- `yout` -- vector/array to write the result to.  If vector it will be
            written like `yout[:]=result`.  If an array the `i`
            argument is also needed and it will be written
            `yout[:,i]=result`.
- `i` -- column of `yout` to write to.

Return: `nothing`
"""
interpolate!(st::AState, tout, yout, i) = hermite_interp!(st, tout, yout, i)
interpolate!(st::AState, tout, yout) = interpolate!(st, tout, yout, 1)

"""
Make dense output using Hermite interpolation of order O(3). Updates
yout in-place.  Only needs y and dy at t-1 and t.

Input

- yout -- inplace y output
- iy   -- column index to write
- tout -- time of requested output
- st::AState -- iteration state variable

Ref: Hairer & Wanner p.190
"""
function hermite_interp!(st::AState, tout, yout, iy)
    y = gety(st)
    dy = getdy(st)
    t = gett(st)

    ypre = gety_pre(st)
    dypre = getdy_pre(st)
    tpre = gett_pre(st)


    if t==tout
        for i=1:length(y)
            yout[i,iy] = y[i]
        end
    elseif tpre==tout
        for i=1:length(y)
            yout[i,iy] = ypre[i]
        end
    else
        dt       = t-tpre
        theta    = (tout-tpre)/dt
        for i=1:length(y)
            yout[i,iy] = ((1-theta)*ypre[i] + theta*y[i] + theta*(theta-1) *
                          ((1-2*theta)*(y[i]-ypre[i]) + (theta-1)*dt*dypre[i] + theta*dt*dy[i]) )
        end
    end
    return nothing
end
