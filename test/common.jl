using DiffEqProblemLibrary, DiffEqBase, ODE

prob = prob_ode_linear
dt=1/2^(4)

algs = [feuler(),rk23(),rk45(),feh78(),ModifiedRosenbrock(),midpoint(),heun(),rk4(),feh45()]

for alg in algs
    sol = solve(prob,alg;dt=dt)
end

prob = prob_ode_2Dlinear

for alg in algs
    if alg != ModifiedRosenbrock() #ODE.jl issues with 2D
      sol = solve(prob,alg;dt=dt)
    end
end

prob = prob_ode_bigfloat2Dlinear

for alg in algs
    if alg != ModifiedRosenbrock() #ODE.jl issues with 2D
        sol = solve(prob,alg;dt=dt)
    end
end
