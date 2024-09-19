using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")

using Plots
using DifferentialEquations

# Create the differential equation system required
function f!(du, u, p::Matrix, t)
    a_11 = p[1,1]
    a_12 = p[1,2]
    a_21 = p[2,1]
    a_22 = p[2,2]

    du[1] = a_11 * u[1] + a_12 * u[2]
    du[2] = a_21 * u[1] + a_22 * u[2]

end

function define_system(A::Matrix, u0::Vector, tspan::Tuple)
    function LinearODE!(du, u, p, t)
        du .= A * u
    end

    return ODEProblem(LinearODE!, u0, tspan)
end

function solve_system(prob::ODEProblem)
    return solve(prob)
end

function plot_solution(sol::ODESolution)
    plot(sol, xlabel = "Time", ylabel = "y", title = "Solution of ODE with A")
end


A = [0.0 1.0; -1.0 0.0]  # System of a simple harmonic oscillator
u0 = [1.0, 0.0]          # Initial conditions
tspan = (0.0, 10.0)      # Time span

prob = define_system(A, u0, tspan)
sol = solve_system(prob)

plot(sol, xlabel = "Time", ylabel = "y", title = "Solution of ODE with A")

# Intial Conditions

#=
u0 = [1.0, 0.0]
A = [0.5 -0.2; 0.1 0.3]
t_span = (0.0, 50.0)

prob = ODEProblem(f!, u0, t_span, A)
sol = solve(prob)

plot(sol, xlabel="Time", ylabel="y", title="Solution of ODE with A")
=#