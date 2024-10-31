using Plots
using DifferentialEquations
using Measurements
using StatsPlots
# Constructing parameters from the data
c = 8 #Number of daily contacts on average
gamma = 1/7 # Daily rate of recovery if it takes 7 days to recover typically
gamma_s = 1/14
p_s = 0.20 # Average probability of severe infection
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month
epsilon = 0.30 # Efficacy of intervention
phi = 0.550 # Proportion of population that will adhere to the intervention

#Using the value of beta that best matches the current data as of 21/10/2024
Beta = measurement(0.035, 0.002)
param_no_int = [c, Beta, gamma, alpha, p_s, gamma_s]
param_int = [c, Beta, gamma, alpha, p_s, gamma_s, epsilon, phi]
# Implement the intervention at 30 days, simulate for another 30
t_span_half_1 = (0, 30)
t_span_half_2 = (30, 120)
pop0_no_int = [S, I, I_s, R]
pop0_int = [S, I, I_s, R]

# Simulate two sets of models, with one set implementing the intervention after day 30

model_no_int1 = ODEProblem(town_SIRS!, pop0_no_int, t_span_half_1, param_no_int)
sol_no_int1 = solve(model_no_int1, saveat = 1)
model_no_int2 = ODEProblem(town_SIRS!, sol_no_int1.u[31], t_span_half_2, param_no_int)
sol_no_int2 = solve(model_no_int2, saveat = 1)

model_int1 = ODEProblem(town_SIRS!, pop0_int, t_span_half_1, param_no_int)
sol_int1 = solve(model_int1, saveat = 1)
model_int2 = ODEProblem(town_SIRS_Intervention!, sol_int1.u[31], t_span_half_2, param_int)
sol_int2 = solve(model_int2, saveat = 1)

# The data of interest is the number of infected, obtain from solution as so
I_model_no_int1 = [u[2].val for u in sol_no_int1.u]
I_model_no_int1_err = [u[2].err for u in sol_no_int1.u]
I_model_no_int2 = [u[2].val for u in sol_no_int2.u]
I_model_no_int2_err = [u[2].err for u in sol_no_int2.u]
append!(I_model_no_int1, I_model_no_int2) 
append!(I_model_no_int1_err, I_model_no_int2_err) 

I_model_int1 = [u[2].val for u in sol_int1.u]
I_model_int1_err = [u[2].err for u in sol_int1.u]
I_model_int2 = [u[2].val for u in sol_int2.u]
I_model_int2_err = [u[2].err for u in sol_int2.u]
append!(I_model_int1, I_model_int2)
append!(I_model_int1_err, I_model_int2_err)

#Append the time together as well
append!(sol_int1.t, sol_int2.t)

#Compare the two models with the data provided up to this point as of 24/10/2024
I_data_d15_d55 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11,
7,
20,
3,
29,
14,
11,
12,
16,
10,
58,
34,
26,
29,
51,
55,
155,
53,
67,
98,
130,
189,
92,
192,
145,
128,
68,
74,
126,
265,
154,
207,
299,273,190,152,276,408,267,462,352]
main = plot()
plot!(main, sol_int1.t, I_model_no_int1, ribbon = I_model_no_int1_err, label = "I Model: No Intervention", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "Intervention at Day 30 vs No Intervention")
plot!(main, sol_int1.t, I_model_int1, ribbon = I_model_int1_err,  label = "I Model: With Intervention at day 30, probability of use = $phi")
plot!(main, range(0, stop = 55, step = 1), I_data_d15_d55, label = "I Data",  seriestype=:scatter)
display(main)