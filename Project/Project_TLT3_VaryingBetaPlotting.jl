#Instead of plotting the errors for each value of Beta, plot the variability in
# the models
using Plots
using DifferentialEquations

#Town Population
N = 6000
I = 1
S = N - I
I_s = 0
R = 0

# Constructing parameters from the data
c = 8 #Number of daily contacts on average
gamma = 1/7 # Daily rate of recovery if it takes 7 days to recover typically
gamma_s = 1/14
p_s = 0.20 # Average probability of severe infection
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month

#Data is operated on a daily basis for 30 days
#We are not provided the first 15 days however
I_data_d15_d30 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,7,20,3,29,14,11,12,16,10,58, 34, 26, 29, 51, 55]
Is_data_d21_d30 = [0, 0, 1, 2, 5,5,5,2,9,4]

t_span = (0, 365)
pop0 = [S, I, I_s, R]
#Beta appears to be between 0.03 and 0.04, find the value between them that minimises
#error
Betas = range(0.032, step = 0.0005, stop = 0.037)

I_models = Array{Any}(undef, length(Betas), t_span[2] + 1)
Is_models = Array{Any}(undef, length(Betas), t_span[2] + 1)
t_sol = 0
for i in 1:length(Betas)
    local param = [c, Betas[i], gamma, alpha, p_s, gamma_s]
    # Create the ODE model of the town
    local model =  ODEProblem(town_SIRS!, pop0, t_span, param)
    local sol = solve(model, saveat = 1)
    # The data of interest is the number of infected, obtain from solution as so
    local I_model = [u[2] for u in sol.u]
    global t_sol = sol.t
    #print(size(I_model))
    I_models[i,:] = I_model
end
# Obtain severe infection model as well 
for i in 1:length(Betas)
    local param = [c, Betas[i], gamma, alpha, p_s, gamma_s]
    # Create the ODE model of the town
    local model =  ODEProblem(town_SIRS!, pop0, t_span, param)
    local sol = solve(model, saveat = 1)
    # The data of interest is the number of infected, obtain from solution as so
    local Is_model = [u[3] for u in sol.u]
    #print(size(I_model))
    Is_models[i,:] = Is_model
end

main = plot(xlabel = "Time(Days)", ylabel = "Number of people in Classed as Infected", title = "SIRS Infected Model with Varying Beta Values")
print(size(I_models))
plot!(main, t_sol, I_data_d15_d30, seriestype=:scatter,  label = "I_data")
for i in 1:length(Betas)
    plot!(main, I_models[i, :], label = "Beta = $(Betas[i])" )
end

display(main)

severe = plot(xlabel = "Time(Days)", ylabel = "Number of people in Classed as Severely Ill", title = "SIRS Severely Ill Model with Varying Beta Values")
print(size(Is_models))
for i in 1:length(Betas)
    plot!(severe, Is_models[i, :], label = "Beta = $(Betas[i])" )
end
display(severe)