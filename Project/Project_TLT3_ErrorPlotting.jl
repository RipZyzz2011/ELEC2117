using Plots
using DifferentialEquations
using Measurements
using Pkg
Pkg.add(path="C:/Users/hamis/.julia/dev/homogenous_SIR_model/")
using homogenous_SIR_model

#This program utilises the Linear Least Squares Estimator (LLSE) in order to find a 
# value for Beta that best correlates with the data

#Not importing correctly for some reason
function error_squares(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
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

#Data is operated on a daily basis for 30 days
#First 15 days of infection data is unknown
I_data_d15_d30 = [11,7,20,3,29,14,11,12,16,10,58,34,26,29,51,55]
Is_data_d21_d30 = [0, 0, 1, 2, 5,5,5,2,9,4]

t_span = (0, 30)
pop0 = [S, I, I_s, R]
#Beta appears to be between 0.02 and 0.04, find the value between them that minimises
#error
Betas = range(0.02, step = 0.00001, stop = 0.1)
#List to store the error sums for each value of beta
b_errors = []
for beta in Betas
    local param = [c, beta, gamma, alpha, p_s, gamma_s]
    # Create the ODE model of the town
    local model =  ODEProblem(town_SIRS!, pop0, t_span, param)
    local sol = solve(model, saveat = 1)
    local I_model = [u[2] for u in sol.u]
    # Compute the LLS error of each model set with the data
    append!(b_errors,error_squares(I_model[15:30], I_data_d15_d30))
end
# Obtain the index with the minimum error
error_index_min = argmin(b_errors)
println("Beta value that gives the smallest error: $(Betas[error_index_min])")
plot(Betas, b_errors, seriestype=:scatter, xlabel = "Beta Values", ylabel = "Least Squares Error", title = "LLSE Error values between Data and Model")
