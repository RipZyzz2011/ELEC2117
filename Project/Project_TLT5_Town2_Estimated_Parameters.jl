using Plots
using DifferentialEquations
using Measurements

#Model the virus in the second town using the LLSE approximated values
town2_Infected_d27_d80 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,21,29, 25, 30, 28, 34, 28, 54, 57, 92, 73, 80, 109, 102, 128, 135, 163, 150, 211, 196, 233, 247, 283, 286, 332, 371, 390, 404, 467, 529, 598, 641, 704, 702, 788, 856, 854, 955, 995, 1065, 1106, 1159, 1217, 1269, 1298, 1328, 1339, 1383, 1431, 1422, 1414, 1485, 1464, 1480]
town2_Severe_d27_d80 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3, 3, 4, 7, 3, 8, 7, 5, 9, 13, 15, 3, 20, 13, 11, 20, 16, 11, 15, 18, 27, 24, 28, 36, 41, 35, 41, 55, 63, 66, 72, 80, 90, 104, 109, 115, 127, 135, 147, 162, 163, 186, 194, 200, 216, 223, 241, 249, 258, 275, 277, 299, 302, 300]

#Town Population
N = 10000
I = 1
S = N - I
I_s = 0
R = 0

c = 8 #Number of daily contacts on average
gamma = 1/7 # Daily rate of recovery if it takes 7 days to recover typically
gamma_s = 1/14
p_s = measurement(0.2, 0.01) # Average probability of severe infection
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month
epsilon = 0.22 # Efficacy of intervention
phi = 0.55 # Proportion of population that will adhere to the intervention

Beta = measurement(0.0371, 0.00)
param_no_int = [c, Beta, gamma, alpha, p_s, gamma_s]
param_int = [c, Beta, gamma, alpha, p_s, gamma_s, epsilon, phi]

# Implement the intervention at day 36, simulate for another 50
t_span_half_1 = (0, 36)
t_span_half_2 = (36, 100)
pop0_no_int = [S, I, I_s, R]
pop0_int = [S, I, I_s, R]

# Simulate two sets of models, with one set implementing the intervention after day 36

model_no_int1 = ODEProblem(town_SIRS!, pop0_no_int, t_span_half_1, param_no_int)
sol_no_int1 = solve(model_no_int1, saveat = 1)
model_no_int2 = ODEProblem(town_SIRS!, sol_no_int1.u[37], t_span_half_2, param_no_int)
sol_no_int2 = solve(model_no_int2, saveat = 1)

model_int1 = ODEProblem(town_SIRS!, pop0_int, t_span_half_1, param_no_int)
sol_int1 = solve(model_int1, saveat = 1)
model_int2 = ODEProblem(town_SIRS_Intervention!, sol_int1.u[37], t_span_half_2, param_int)
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

# Obtain the modelling data of the severe illness
Is_model_no_int1 = [u[3].val for u in sol_no_int1.u]
Is_model_no_int1_err = [u[3].err for u in sol_no_int1.u]
Is_model_no_int2 = [u[3].val for u in sol_no_int2.u]
Is_model_no_int2_err = [u[3].err for u in sol_no_int2.u]
append!(Is_model_no_int1, Is_model_no_int2) 
append!(Is_model_no_int1_err, Is_model_no_int2_err)

Is_model_int1 = [u[3].val for u in sol_int1.u]
Is_model_int2 = [u[3].val for u in sol_int2.u]
Is_model_int1_err = [u[3].err for u in sol_int1.u]
Is_model_int2_err = [u[3].err for u in sol_int2.u]
append!(Is_model_int1, Is_model_int2)
append!(Is_model_int1_err, Is_model_int2_err)

#Append the time together as well
append!(sol_int1.t, sol_int2.t)

#Plot the data of town 2 alongside the model 
main = plot()
plot!(main, sol_int1.t, I_model_no_int1, ribbon = I_model_no_int1_err, label = "I Model: No Intervention", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "Infected Popultation of Town 2")
plot!(main, sol_int1.t, I_model_int1, ribbon = I_model_int1_err,  label = "I Model: With Intervention at day 36, probability of use = $phi")
#plot!(main, range(0, stop = 55, step = 1), I_data_d15_d55, label = "I Data",  seriestype=:scatter)
plot!(main, range(0, stop = 80, step = 1), town2_Infected_d27_d80, label = "Infected population in town 2 data as of 31/10",  seriestype=:scatter)
display(main)


severe = plot()
plot!(severe, sol_int1.t, Is_model_no_int1, ribbon = Is_model_no_int1_err, label = "Is Model: No Intervention", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "Severe Illness Population of Town 2")
plot!(severe, sol_int1.t, Is_model_int1, ribbon = Is_model_int1_err,  label = "Is Model: With Intervention at day 36, probability of use = $phi")
plot!(severe, range(0, stop = 80, step = 1),town2_Severe_d27_d80, label = "Severe Illness Data in Town 2 as of 31/10",  seriestype=:scatter)
display(severe)
display(main)