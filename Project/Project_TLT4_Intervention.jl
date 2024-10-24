# Recalibrating our model to include the intervention method and its
# associated rate

#SIR model that now incorporates re-infection and severe illness state, as well
# as an intervention rate and probability

#=
  It is anticipated that this intervention will reduce the probability of transmission from an infected person
   to a susceptible person by about 30% (this is the efficacy of the intervention). 
   The Department of Health expects approximately 
   80% of the population will comply with their directive to implement the public health intervention.
=#
function town_SIRS_Intervention!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta_c, gamma, alpha, p_s, gamma_s, epsilon, phi = param
    #gamma: Probability of recovering from infection each day
    S, I, Is, R = pop
    #Lambda: Force of Infection
    lambda = c * (1 - epsilon * phi) * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S + alpha * R# dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    #Severe infection
    dpop[3] = gamma * p_s * I - gamma_s * Is 
    dpop[4] = (1 - p_s) * gamma * I + gamma_s * Is - alpha * R # dR = gamma * R

end

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
p_s = 0.20 # Average probability of severe infection
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month
epsilon = 0.30 # Efficacy of intervention
phi = 0.80 # Proportion of population that will adhere to the intervention

#Using the value of beta that best matches the current data as of 21/10/2024
Beta = 0.03513
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
I_model_no_int1 = [u[2] for u in sol_no_int1.u]
I_model_no_int2 = [u[2] for u in sol_no_int2.u]
append!(I_model_no_int1, I_model_no_int2) 

I_model_int1 = [u[2] for u in sol_int1.u]
I_model_int2 = [u[2] for u in sol_int2.u]
append!(I_model_int1, I_model_int2)

#Append the time together as well
append!(sol_int1.t, sol_int2.t)
main = plot()
plot!(main, sol_int1.t, I_model_no_int1, label = "I Model: No Intervention", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "Intervention at Day 30 vs No Intervention")
plot!(main, sol_int1.t, I_model_int1,  label = "I Model: With Intervention at day 30")

display(main)