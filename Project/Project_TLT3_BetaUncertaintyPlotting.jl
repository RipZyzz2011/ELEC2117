#Same idea as the other TLT 3 file, however instead of plotting a range
# beta is now going to have an uncertainty attached
using Plots
using DifferentialEquations
using Measurements
using Pkg
Pkg.add(["Measurements", "StatsPlots"])
using homogenous_SIR_model

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
I_data_d15_d30 = [11,7,20,3,29,14,11,12,16,10,58, 34, 26, 29, 51, 55]
Is_data_d21_d30 = [ 1, 2, 5,5,5,2,9,4]

t_span = (0, 30)
pop0 = [S, I, I_s, R]

# Set beta to 0.035 with an uncertainty of +/- 0.002
beta = measurement(0.035, 0.003)
R_0 = c * beta / gamma
param = [c, beta, gamma, alpha, p_s, gamma_s]
# Create the ODE model of the town
model =  ODEProblem(town_SIRS!, pop0, t_span, param)
sol = solve(model, saveat = 1)
# The data of interest is the number of infected, obtain from solution as so
I_model_mean = [(u[2].val) for u in sol.u]
I_model_err = [(u[2].err) for u in sol.u]
Is_model_mean = [(u[3].val) for u in sol.u]
Is_model_err = [(u[3].err) for u in sol.u]



println("Beta is approximately $beta from the data, and R_0 is approximately $R_0")
plot(sol.t, I_model_mean, ribbon = I_model_err, label = "I_Model", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIRS Infected Model vs Infected Population Data")
plot!(range(15, step = 1, stop = 30), I_data_d15_d30, seriestype=:scatter, label = "I_Data") 
plot(sol.t, Is_model_mean, ribbon = Is_model_err, label = "Is_Model", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIRS Severe Illness Model vs Severe Illness Data")
plot!(range(21, step = 1, stop = 30), Is_data_d21_d30, seriestype=:scatter, label = "Is_Data")     

