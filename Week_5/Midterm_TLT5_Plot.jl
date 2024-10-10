using Pkg
using Plots
Pkg.add(path="C:/Users/hamis/.julia/dev/homogenous_SIR_model/")
using homogenous_SIR_model

#TLT Level 5, now 15% of the population is already immune
#Same implementation as TLT4, however now 15% of the population starts as recovered
#In this implementation, recovered people cannot get sick again
#From the problem description, define the parameters
N = 10000 #Town population
#Initial town conditions
S = 0.845 * N
I = 0.0005 * N
R = 0.15 * N

pop0 = [S, I, R] #Compile population into an array

c = 15 # Number of daily contacts
#From use of the error mapping, 0.044 seems to be the best value of Beta
Beta = 0.061 # Rate of infection per contact
gamma = 0.10 # Daily recovery rate

param = [c, Beta, gamma]

#System operates over a 30-day period
tspan = (0, 30)
# Generate the model using the homogenous_SIR_model package, using the herd immunity
# criteria model given the questions request
model = define_town_model(:foi, param, pop0, tspan)
# solve the system to allow for plotting
sol = solve_system(model)
# The data of interest is the number of infected, obtain from solution as so
I_model = [u[2] for u in sol.u]

# Additionally, based on the parameters, obtain the herd immunity threshold
# Obtain using the reproduction number based on the parameters
R_0 = c * 1/gamma * Beta
#Herd immunity threshold:
p_c = 1 - 1/R_0

#Plot the model given the new immunity condition
#u[1] = susceptible, u[2] = infected, u[3] = recovered
#plot_model_solution(sol)

# plot against the given data to compare model
plot(sol.t, I_model, label = "I_Model", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "TLT5: SIR Model vs Data")
plot!(town_data, seriestype=:scatter, label = "I_Data")   
