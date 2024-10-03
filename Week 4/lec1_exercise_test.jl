using Test
include("lec1_exercise.jl")

@testset "townSIR_two_param_set" begin
    #Test that no one gets infected when a population has I = 0
    S = 100
    I = 0
    R = 0
    N_initial = S + I + R
    lambda = 0.1
    gamma = 0.1

    t_span = (0, 50)

    model = define_town_model((lambda, gamma), [S, I, R], t_span)
    solve_system(model)
    
    @test I == 0

    #Test that the population remains the same size that it begun at
    S = 90
    I = 10
    R = 0
    N_initial = S + I + R
    lambda = 0.1
    gamma = 0.1

    t_span = (0, 50)

    model = define_town_model((lambda, gamma), [S, I, R], t_span)
    solve_system(model)
    
    @test (S + I + R) == N_initial
end