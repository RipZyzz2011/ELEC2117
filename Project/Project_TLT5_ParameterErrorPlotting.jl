# Initially assume that virus behaves the same as the previous town
# however the intervention has varied effectiveness 
# The purpose of this script is to vary epsilon (intervention efficacy) and plot 
# its error in relation to the data
using Plots
using DifferentialEquations
using Measurements
#using homogenous_SIR_model

#Calculate the sum of the error at each datapoint between the real values and the model values
#Squares error between values, useful for evaluating the efficacy of the beta parameter
function error_squares(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end

# Investigate and model the pathogen in the second town.
# First detected on day 27 relative to town 1, intervention begun on day 36 relative to town 1
town2_Infected_d27_d80 = [21,29, 25, 30, 28, 34, 28, 54, 57, 92, 73, 80, 109, 102, 128, 135, 163, 150, 211, 196, 233, 247, 283, 286, 332, 371, 390, 404, 467, 529, 598, 641, 704, 702, 788, 856, 854, 955, 995, 1065, 1106, 1159, 1217, 1269, 1298, 1328, 1339, 1383, 1431, 1422, 1414, 1485, 1464, 1480]
town2_Severe_d27_d80 = [3, 3, 4, 7, 3, 8, 7, 5, 9, 13, 15, 3, 20, 13, 11, 20, 16, 11, 15, 18, 27, 24, 28, 36, 41, 35, 41, 55, 63, 66, 72, 80, 90, 104, 109, 115, 127, 135, 147, 162, 163, 186, 194, 200, 216, 223, 241, 249, 258, 275, 277, 299, 302, 300]

#Initially begin modelling the virus with the parameters of the first town's virus
#Town Population
N = 10000
I = 1
S = N - I
I_s = 0
R = 0

c = 8 #Number of daily contacts on average
gamma = 1/7 # Daily rate of recovery if it takes 7 days to recover typically
gamma_s = 1/14
p_s = measurement(0.20, 0.05) # Average probability of severe infection
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month
phi = 0.55 # Proportion of population that will adhere to the intervention

Beta = measurement(0.035, 0.002)

#Vary the efficacy of the intervention whilst maintaining every other parameter
Epsilons = range(start = 0, step = 0.01, stop = 1)
t_span = (0, 80)
pop0 = [S, I, I_s, R]

e_errors_infected = []
e_errors_severe = []

for epsilon in Epsilons
    local param_int = [c, Beta, gamma, alpha, p_s, gamma_s, epsilon, phi]
     # Create the ODE model of the town
     local model =  ODEProblem(town_SIRS_Intervention!, pop0, t_span, param_int)
     local sol = solve(model, saveat = 1)
     # The data of interest is the number of infected, obtain from solution as so
     local I_model = [u[2].val for u in sol.u]
     local Is_model = [u[3].val for u in sol.u]
     append!(e_errors_infected, error_squares(I_model[27:81],town2_Infected_d27_d80))
     append!(e_errors_severe, error_squares(Is_model[27:81], town2_Severe_d27_d80))

end
# Obtain the indices with the minimum error
error_infected_index_min = argmin(e_errors_infected)
error_severe_index_min = argmin(e_errors_severe)
println("Epsilon value that gives the smallest infected error: $(Epsilons[error_infected_index_min])")
println("Epsilon value that gives the smallest severe illness error: $(Epsilons[error_severe_index_min])")
errors = plot()
plot!(errors, Epsilons, e_errors_infected, seriestype=:scatter, xlabel = "Epsilon Values", ylabel = "Least Squares Error", label = "Infected People Errors")
plot!(errors, Epsilons, e_errors_severe, seriestype=:scatter, label = "Severe Illness Errors")

######
######
######
# Repeat the process but now disregard intervention, and simply vary the transmission
# parameters
                                                                        #^
#Use the value of epsilon + uncertainty determined from the test above   |
epsilon = measurement(0.15, 0.03)
#Vary beta between the range of reasonable sensitivity
Betas = range(0.02, step = 0.00001, stop = 0.04)
b_errors_infected = []
b_errors_severe = []
for beta in Betas
    local param_int = [c, beta, gamma, alpha, p_s, gamma_s]
     # Create the ODE model of the town
     local model =  ODEProblem(town_SIRS!, pop0, t_span, param_int)
     local sol = solve(model, saveat = 1)
     # The data of interest is the number of infected, obtain from solution as so
     local I_model = [u[2].val for u in sol.u]
     local Is_model = [u[3].val for u in sol.u]
     append!(b_errors_infected, error_squares(I_model[27:81],town2_Infected_d27_d80))
     append!(b_errors_severe, error_squares(Is_model[27:81], town2_Severe_d27_d80))

end
# Obtain the indices with the minimum error
beta_infected_min = argmin(b_errors_infected)
beta_severe_min = argmin(b_errors_severe)
println("Beta value that gives the smallest infected error: $(Betas[beta_infected_min])")
println("Beta value that gives the smallest severe illness error: $(Betas[beta_severe_min])")
beta_plotting = plot()
plot!(beta_plotting, Betas, b_errors_infected, seriestype=:scatter, xlabel = "Beta Values", ylabel = "Least Squares Error", label = "Infected People Errors")
plot!(beta_plotting, Betas, b_errors_severe, seriestype=:scatter, label = "Severe Illness Errors")
display(beta_plotting)
