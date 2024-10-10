using Plots
using DifferentialEquations
#using homogenous_SIR_model

#Calculate the sum of the error at each datapoint between the real values and the model values
#Squares error between values, useful for evaluating the efficacy of the beta parameter
function error_squares(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end
#SIR model that now incorporates re-infection and severe illness state
function town_SIRS!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta_c, gamma, alpha, p_s, gamma_s = param
    #gamma: Probability of recovering from infection each day
    S, I, Is, R = pop
    lambda = c * Beta_c * I / N
    R_0 = c * 1/gamma * Beta_c # Reproduction number
    
    dpop[1] = -lambda * S + alpha * R# dS = -lambda*S
    dpop[2] = lambda * S - gamma * I # dI = lambda * S - gamma * R
    #Severe infection
    dpop[3] = gamma * p_s * I - gamma_s * Is 
    dpop[4] = (1 - p_s) * gamma * I + gamma_s * Is - alpha * R # dR = gamma * R

end

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

#Data is operated on a daily basis for 25 days
I_data_d14_d25 = [11,7,20,3,29,14,11,12,16,10,58]
Is_data_d21_d25 = [0, 0, 1, 2, 5]

t_span = (0, 25)
pop0 = [S, I, I_s, R]
#Beta appears to be between 0.03 and 0.04, find the value between them that minimises
#error
Betas = range(0.02, step = 0.00001, stop = 0.04)
b_errors = []
for beta in Betas
    local param = [c, beta, gamma, alpha, p_s, gamma_s]
    # Create the ODE model of the town
    local model =  ODEProblem(town_SIRS!, pop0, t_span, param)
    local sol = solve(model, saveat = 1)
    # The data of interest is the number of infected, obtain from solution as so
    local I_model = [u[2] for u in sol.u]
    append!(b_errors,error_squares(I_model[14:26], I_data_d14_d25))
end
# Obtain the index with the minimum error
error_index_min = argmin(b_errors)
println("Beta value that gives the smallest error: $(Betas[error_index_min])")
plot(Betas, b_errors, seriestype=:scatter, xlabel = "Beta Values", ylabel = "Least Squares Error")
