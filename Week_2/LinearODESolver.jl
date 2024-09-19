module LinearODESolver
using Plots
using DifferentialEquations
using Test

export define_system, solve_system, plot_solution

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

@testset "LinearODESolver Tests" begin

    A = [0.0 1.0; -1.0 0.0]  # System of a simple harmonic oscillator
    u0 = [1.0, 0.0]          # Initial conditions
    tspan = (0.0, 10.0)      # Time span

    prob = define_system(A, u0, tspan)
    sol = solve_system(prob)

    println(typeof(sol))
    
    @test typeof(sol) == ODESolution
    @test isapprox(sol[1,13], cos(10), atol=1e-2)  # Verify result at t = 10.0
    @test isapprox(sol[2,13], -sin(10), atol=1e-2) # Check solution of oscillator

    # 3. Test plot_solution (just test that it runs without error)
    @test plot_solution(sol) isa Plots.Plot

end

end


