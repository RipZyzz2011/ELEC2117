using Plots
using DifferentialEquations

#Town size
#lambda = c(t) * Beta_c(t) * I(t) / N
# c: Number of daily contacts
# Beta_c: Chance of contracting disease from infected contact
# I(t) / N: Chance of infectious contact given homogenous population
function town_SIR!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta_c, gamma, p_c = param
    #gamma: Probability of recovering from infection each day
    S, I, R = pop
    lambda = c * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S # dS = -lambda*S
    dpop[2] = gamma * R_0 * S * I / N - gamma * I # dI = lambda * S - gamma * R
    dpop[3] = gamma * I # dR = gamma * R

    p_c = 1 - 1/R_0 # Herd immunity threshold

end

# Run the model with some initial conditions
pop0 = [4999, 1, 0]
tspan = (0.0, 22)
c = 10
Beta_c = 0.03
gamma = 0.1
p_c = 0
param = [c, Beta_c, gamma, p_c] #c, Beta_c, gamma
prob = ODEProblem(town_SIR!, pop0, tspan, param)

sol = solve(prob)
data = [1.0, 0.0, 5.0, 12.0, 0.0, 12.0, 0.0, 12.0, 11.0, 13.0, 0.0, 17.0, 41.0, 27.0, 20.0, 41.0, 47.0, 61.0, 76.0, 113.0, 158.0]

plot(sol, xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIR Model")
   