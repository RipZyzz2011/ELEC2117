using Pkg
Pkg.add(["Measurements", "StatsPlots", "Plots", "DifferentialEquations"])
#Change this path to the relevant local Github repositry
Pkg.add(path = "C:/Users/hamis/.julia/dev/homogenous_SIR_model/")
using Plots
using DifferentialEquations
using homogenous_SIR_model

# This script is the starting point in accurately modelling the outbreak in the town
# as such it forgoes incorportating uncertainties and serves to only answer the questions
# asked in TLT 2

#Model the progression of the disease for the first 30 days of the outbreak
# alongside the data provided
#Town Population, assume that there was a single person infected on day 1
N = 6000
I = 1
S = N - I
I_s = 0
R = 0


# Constructing parameters from the data
c = 8 #Number of daily contacts on average
gamma = 1/7 # Daily rate of recovery if it takes 7 days to recover typically
gamma_s = 1/14 # Daily rate of recovery from severe illness
p_s = 0.20 # Average probability of severe infection
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month

#Data is operated on a daily basis for 25 days
I_data_d14_d25 = [11,7,20,3,29,14,11,12,16,10,58]
Is_data_d21_d25 = [0, 0, 1, 2, 5]
#Using the least-squares error method, the beta that provides the smallest error
# with the current data is 0.03697
Beta = 0.03697
#Estimate R0, not accounting for the 
R_0 = c * 1/gamma * Beta

param = [c, Beta, gamma, alpha, p_s, gamma_s]
t_span = (0, 25)
pop0 = [S, I, I_s, R]

# Create the ODE model of the town
model =  ODEProblem(town_SIRS!, pop0, t_span, param)
sol = solve(model, saveat = 1)
# The data of interest is the number of infected, obtain from solution as so
I_model = [u[2] for u in sol.u]

println("Beta is approximately $Beta from the data, and R_0 is approximately $R_0")
plot(sol.t, I_model, label = "I_Model", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIRS Model vs Data")
plot!(range(14, step = 1, stop = 25),I_data_d14_d25, seriestype=:scatter, label = "I_Data")   


