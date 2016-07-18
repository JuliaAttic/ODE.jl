# unpack macro from Parameters.jl
#########
"""
Unpacks fields from any datatype (no need to create it with @with_kw):

```julia
type A
    a
    b
end
aa = A(3,4)
@unpack aa: a,b
# is equivalent to
a = aa.a
b = aa.b
```
"""
macro unpack(arg)
    v, up = parse_pack_unpack(arg)
    out = quote end
    for u in up
        push!(out.args, :($u = $v.$u))
    end
    return esc(out)
end
## @pack and @unpack are independent of the datatype
function parse_pack_unpack(arg)
    if isa(arg, Symbol); error("Need format `t: a, ...`") end
    h = arg.head
    if !(h==:(:) || h==:tuple); error("Need format `t: a, ...`") end
    # var-name of structure
    v = h==:tuple ? arg.args[1].args[1] : arg.args[1]
    # vars to unpack
    up = Any[]
    if h==:tuple
        append!(up, arg.args[2:end])
        push!(up, arg.args[1].args[2])
    else
        push!(up, arg.args[2])
    end
    return v, up # variable holding the structure, variables to (un)-pack
end


##########
# Conversion to in-place function
##########
inplace(F) = (t,y,dy) -> copy!(dy,F(t,y))
outofplace(F!) = (t,y) -> (dy=similar(y); F!(t,y,dy); dy)

########
# Jacobians
#######
# from @pwl's PR49

# generate a jacobian using ForwardDiff
"Automatic differentiation Jacobian"
ad_jacobian(F) = (t,y)->ForwardDiff.jacobian(y->F(t,y),y)
"Finite difference Jacobian"
fd_jacobian(F) = (t,y)->_fd_jacobian(F,y,t)
function _fd_jacobian(F, y, t)
    fty = F(t, y)
    ly = max(length(y),1)
    dFdy = zeros(eltype(y), ly, ly)
    for j = 1:ly
        # The 100 below is heuristic
        dy = zeros(eltype(y), ly)
        dy[j] = (y[j] .+ (y[j]==0))./100
        dFdy[:,j] = (F(t,y+dy)-fty)./dy[j]
    end
    return dFdy
end

###
# Initial step size guess
###
# "Solving Ordinary Differential Equations I" by Hairer et al., p.169
function hinit(F, x0, t0, tend, p, reltol, abstol)
    # Returns first step, direction of integration and F evaluated at t0
    tdir = sign(tend-t0)
    tdir==0 && error("Zero time span")
    tau = max(reltol*norm(x0, Inf), abstol)
    d0 = norm(x0, Inf)/tau
    f0 = F(t0, x0)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    # perform Euler step
    x1 = x0 + tdir*h0*f0
    f1 = F(t0 + tdir*h0, x1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 1e-15
        h1 = max(1e-6, 1e-3*h0)
    else
        pow = -(2. + log10(max(d1, d2)))/(p + 1.)
        h1 = 10.^pow
    end
    return tdir*min(100*h0, h1, tdir*(tend-t0)), tdir, f0
end

####
# Misc
##
tdir(ode) = sign(ode.tspan[2]-ode.tspan[1])


# function dict2type(di::Dict{Symbol}, T::DataType)
#     di = Dict{Symbol,Any}()
#     for n in fieldnames(dt)
#         di[n] = getfield(dt, n)
#     end
#     di
# end
