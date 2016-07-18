# https://github.com/JuliaLang/ODE.jl/issues/104
# Works

if !isdefined(:IVP)
    include("ODE.jl")
end
function f(t, y)
    # Extract the components of the y vector
    (x, v) = y

    # Our system of differential equations
    x_prime = v
    v_prime = -x

    # Return the derivatives as a vector
    [x_prime; v_prime]
end

y0 = [0.0; 0.1]
tsteps = 0.:-.1:-4pi;

@time tout,yout = ode(RK4, f, y0, tsteps)
@time tout,yout = ode(Rosenbrock23s{:v2}, f, y0, tsteps)
