using Plots
using DifferentialEquations

#Town size
#lambda = c(t) * Beta_c(t) * I(t) / N
# c: Number of daily contacts
# Beta_c: Chance of contracting disease from infected contact
# I(t) / N: Chance of infectious contact given homogenous population
function town_SIR!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta_c, gamma = param
    #gamma: Probability of recovering from infection each day
    S, I, R = pop
    lambda = c * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S # dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    dpop[3] = gamma * I # dR = gamma * R

end

# Run the model with some initial conditions
pop0 = [99, 1, 0]
tspan = (0.0, 50.0)
c = 10
Beta_c = 0.03
gamma = 0.05
param = [c, Beta_c, gamma] #c, Beta_c, gamma
prob = ODEProblem(town_SIR!, pop0, tspan, param)

sol = solve(prob)

plot(sol, xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIR Model")    