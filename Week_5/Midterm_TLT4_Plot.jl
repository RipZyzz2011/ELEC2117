# The model is packaged
using Pkg
using Plots
Pkg.add(path="C:/Users/hamis/.julia/dev/homogenous_SIR_model/")
using homogenous_SIR_model

# TLT 4: Comparing the model to the actual data of the outbreak in the town
#30 day period data
town_data = [5, 10, 19, 37, 71, 136,
260, 486, 882, 1516, 2399, 3407, 4300, 4882, 5116, 5080, 4875, 4582, 4251, 3913, 3583, 3271,
2979, 2708, 2460, 2233, 2026, 1837, 1665, 1509]

#= Parameter description
The Department of Health has reported the early stages of an outbreak of a contagious disease
in a town with a population of 10,000 people. The lab is asked to analyse the potential spread
of the disease. Initial data shows that 0.05% of the population is currently infected, while the
rest are susceptible. Based on previous outbreaks of similar diseases, the Department estimates
that each infected individual is likely to infect 3% of the susceptible individuals they come
into contact with. Additionally, it is assumed that each person has an average of 15 contacts
per day. Furthermore, 10% of infected individuals are expected to recover daily.
=#

#From the problem description, define the parameters
N = 10000 #Town population
#Initial town conditions
S = 0.995 * N
I = 0.0005 * N
R = 0

pop0 = [S, I, R] #Compile population into an array

c = 15 # Number of daily contacts
#From use of the error mapping, 0.037 seems to be the best value of Beta
Beta = 0.05 # Rate of infection per contact
gamma = 0.10 # Daily recovery rate

param = [c, Beta, gamma]

#System operates over a 30-day period
tspan = (0, 30)
# Generate the model using the homogenous_SIR_model package, using the herd immunity
# criteria model given the questions request
model = define_town_model(:herd, param, pop0, tspan)
# solve the system to allow for plotting
sol = solve_system(model)
# The data of interest is the number of infected, obtain from solution as so
I_model = [u[2] for u in sol.u]

# Additionally, based on the parameters, obtain the herd immunity threshold
# Obtain using the reproduction number based on the parameters
R_0 = c * 1/gamma * Beta
#Herd immunity threshold for this given population:
p_c = (1 - 1/R_0) * N
println("Number of people that must be infected for herd immunity: $p_c")
# plot against the given data to compare model
#plot_model_solution(sol)
plot(sol.t, I_model, label = "I_Model", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "TLT4: SIR Model vs Data")
plot!(town_data, seriestype=:scatter, label = "I_Data")   

