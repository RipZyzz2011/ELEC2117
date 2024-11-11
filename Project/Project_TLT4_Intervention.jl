using Plots
using DifferentialEquations
using Measurements
using Pkg
Pkg.add(path="C:/Users/hamis/.julia/dev/homogenous_SIR_model/")
using homogenous_SIR_model
# Recalibrating our model to include the intervention method and its
# associated rate

#SIR model that now incorporates re-infection and severe illness state, as well
# as an intervention rate and probability

# Compare the impact of intervention after day 30 with no intervention
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
p_s = measurement(0.20, 0.05) # Average probability of severe infection
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month
epsilon = 0.30 # Efficacy of intervention
phi = 0.55 # Proportion of population that will adhere to the intervention

#Using the value of beta that best matches the current data as of 21/10/2024
Beta = measurement(0.035, 0.002)
param_no_int = [c, Beta, gamma, alpha, p_s, gamma_s]
param_int = [c, Beta, gamma, alpha, p_s, gamma_s, epsilon, phi]
# Implement the intervention at 30 days, simulate for another 30
t_span_half_1 = (0, 30)
t_span_half_2 = (30, 365)
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
main = plot()
plot!(main, sol_int1.t, I_model_no_int1, ribbon = I_model_no_int1_err, label = "I Model: No Intervention", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "Intervention at Day 30 vs No Intervention")
plot!(main, sol_int1.t, I_model_int1, ribbon = I_model_int1_err, label = "I Model: With Intervention at day 30")

display(main)

#Create another plot for severe illness
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

severe = plot()
plot!(severe, sol_int1.t, Is_model_no_int1, label = "Is Model: No Intervention", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "Comparing Severe Illness numbers with Intervention at Day 30 vs No Intervention")
plot!(severe, sol_int1.t, Is_model_int1,  label = "Is Model: With Intervention at day 30")
display(severe)

#Show both on the one plot
plot!(main, sol_int1.t, Is_model_no_int1, ribbon = Is_model_no_int1_err, label = "Is Model: No Intervention")
plot!(main, sol_int1.t, Is_model_int1, ribbon = Is_model_int1_err,  label = "Is Model: With Intervention at day 30")
display(main)