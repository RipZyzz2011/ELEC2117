using Plots
using DifferentialEquations

#Town size
#lambda = c(t) * Beta_c(t) * I(t) / N
# c: Number of daily contacts
# Beta_c: Chance of contracting disease from infected contact
# I(t) / N: Chance of infectious contact given homogenous population
function town_SIR!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta_c, gamma, alpha = param
    #gamma: Probability of recovering from infection each day
    S, I, R = pop
    lambda = c * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S + alpha * R# dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    dpop[3] = gamma * I - alpha * R # dR = gamma * R

    p_c = 1 - 1/R_0 # Herd immunity threshold

end

#SIR model that now incorporates re-infection and severe illness state
function town_SIRS!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta_c, gamma, alpha, p_s = param
    #gamma: Probability of recovering from infection each day
    S, I, IS, R = pop
    lambda = c * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S + alpha * R# dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    #Severe infection
    dpop[3] = gamma * ps * I - sgamma
    dpop[4] = gamma * I - alpha * R # dR = gamma * R
    

    p_c = 1 - 1/R_0 # Herd immunity threshold

end

#Calculate the sum of the error at each datapoint between the real values and the model values
function error(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end

# Run the model with some initial conditions
pop0 = [4999, 1, 0]
tspan = (0.0, 22.0)
c = 10
Beta_c = 0.03519
gamma = 0.1
alpha = 0.1
param = [c, Beta_c, gamma, alpha] #c, Beta_c, gamma
prob = ODEProblem(town_SIR!, pop0, tspan, param)

sol = solve(prob, saveat = 1)
# Actual infection data given to us by Dept' of Health
data = [1.0, 0.0, 5.0, 12.0, 0.0, 12.0, 0.0, 12.0, 11.0, 13.0, 0.0, 17.0, 41.0, 27.0, 20.0, 41.0, 47.0, 61.0, 76.0, 113.0, 158.0]

#Print the error for a given beta values
error_beta = error(I, data)
println("Error for Beta = $Beta_c: $error_beta")


println(length(I))
plot(sol.t, I, label = "model", xlabel = "Time(Days)", ylabel = "Number of people in Category", title = "SIR Model")
plot!(data, seriestype=:scatter, label = "data")   

