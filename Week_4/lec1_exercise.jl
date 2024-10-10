using DifferentialEquations
export define_town_model, solve_system, plot_model_solution

function define_town_model(parameters::Tuple, initial_pop::Vector, t_span::Tuple)
    #=Function: town_SIR
    # Creates a Susceptible, Infected, Recovered (SIR) mathematical model of a town
    on the assumption and basis of a homogenous population, and that recovered people
    aren't susceptible to infection again

    dpop: vector containing the rate of change of each population category
    pop: vector containing the current number of each population category
    param: vector containing the parameters of the differential DifferentialEquations
    t: Beginning and starting time that the model operates from in days

    Each version of the function operates differently based on the size of the param vector, 
    =#

    # Operates with a specified lambda
    function town_SIR!(dpop, pop, parameter::Tuple{2}, t)
        N = sum(pop)
        lambda, gamma = parameter 
    
        S, I, R = pop
    
        dpop[1] = -lambda * S # dS = -lambda*S
        dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
        dpop[3] = gamma * I # dR = gamma * R
    
    end


    # Operates with the force of infection equation for lambda
    function town_SIR!(dpop, pop, param::Tuple{3}, t)
        N = sum(pop)
        c, Beta_c, gamma = param 
        #c: Number of daily contacts
        #Beta_c: Chance of contracting disease from infected contact
        #gamma: Probability of recovering from infection each day
    
        S, I, R = pop
        lambda = c * Beta_c * I / N #I(t) / N: Chance of infectious contact given homogenous population 
        #R_0 = c * 1/gamma * Beta_c # Reproduction number
    
        dpop[1] = -lambda * S # dS = -lambda*S
        dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
        dpop[3] = gamma * I # dR = gamma * R
    
    end

    #= Operates using R0 and the herd immunity threshold criteria
    function town_SIR!(dpop, pop, par::Vector{3}, t)
        N = sum(pop)
        c, Beta_c, gamma = par 
        #c: Number of daily contacts
        #Beta_c: Chance of contracting disease from infected contact
        #gamma: Probability of recovering from infection each day
    
        S, I, R = pop
        lambda = c * Beta_c * I / N #I(t) / N: Chance of infectious contact given homogenous population 
        R_0 = c * 1/gamma * Beta_c # Reproduction number
    
        dpop[1] = -lambda * S # dS = -lambda*S
        dpop[2] = gamma * I * (R0 * S/N - 1)
        dpop[3] = gamma * I # dR = gamma * R
    
    end
    =#
    

    return ODEProblem(town_SIR!, initial_pop, t_span, parameters)
end


#= 
Function: solve_system
Produces the solution to the model created by define_town_model
=#
function solve_system(prob::ODEProblem)
    return solve(prob)
end
#=
Function: plot_model_solution
Plots the solved model over the time span it was solved from
=#

function plot_model_solution(sol::ODESolution)
    plot(sol, xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIR Model")
end