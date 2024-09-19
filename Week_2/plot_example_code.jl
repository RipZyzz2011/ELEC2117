using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")

using Plots
using DifferentialEquations

function f!(du, u, p, t)
    du[1] = -2u[1] + 1
end

u0 = [0.0]
tspan = (0.0, 5.0)
prob = ODEProblem(f!, u0, tspan)

sol = solve(prob)

plot(sol, xlabel = "Time", ylabel = "y", title = "Solution of ODE")