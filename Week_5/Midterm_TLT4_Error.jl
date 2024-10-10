# Find the Beta parameter that best matches the data
# Make use of the error function 
using Pkg
using Plots
Pkg.add(path="C:/Users/hamis/.julia/dev/homogenous_SIR_model/")
using homogenous_SIR_model
#Calculate the sum of the error at each datapoint between the real values and the model values
#Squares error between values, useful for evaluating the efficacy of the beta parameter
function error_squares(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end

# Same parameters and operation as the plot to begin with 
town_data = [5, 10, 19, 37, 71, 136,
260, 486, 882, 1516, 2399, 3407, 4300, 4882, 5116, 5080, 4875, 4582, 4251, 3913, 3583, 3271,
2979, 2708, 2460, 2233, 2026, 1837, 1665, 1509]
#From the problem description, define the parameters
N = 10000 #Town population
#Initial town conditions
S = 0.995 * N
I = 0.005 * N
R = 0

pop0 = [S, I, R] #Compile population into an array

c = 15 # Number of daily contacts
# Span beta over a set surrounding the provided value
Betas = [0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035,
0.036, 0.037, 0.038, 0.039, 0.04]
gamma = 0.10 # Daily recovery rate

error_list = []

#System operates over a 30-day period
tspan = (0, 30)
# Evaluate the error for each value of beta
for beta in Betas
    local param = [c, beta, gamma]
    # Generate the model using the homogenous_SIR_model package, using the herd immunity
    # criteria model given the questions request
    local model = define_town_model(:herd, param, pop0, tspan)
    # solve the system to allow for plotting
    local sol = solve_system(model)
    # The data of interest is the number of infected, obtain from solution as so
    local I = [u[2] for u in sol.u]

    b_error = error_squares(I, town_data)
    append!(error_list, b_error)
end
print(error_list)

