using Test
using LinearODESolver

using DifferentialEquations

# Test the LinearODESolver module
@testset "LinearODESolver Tests" begin

    # Define test system
    A = [0.0 1.0; -2.0 -3.0]
    u0 = [1.0, 0.0]
    tspan = (0.0, 10.0)

    # 1. Test define_system
    prob = define_system(A, u0, tspan)
    @test typeof(prob) == ODEProblem

    # 2. Test solve_system
    sol = solve_system(prob)
    @test typeof(sol) == ODESolution

    # Check solution at a specific point (e.g., t=10)
    known_solution = [2.718; -1.234]  # Example: Add known values if available
    @test isapprox(sol(10), known_solution, atol=0.01)

    # 3. Test plot_solution (just test that it runs without error)
    @test plot_solution(sol) isa Plots.Plot

end