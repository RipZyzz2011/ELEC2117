# Purpose of this file is to perform LLSE whilst varying several parameters at once to find
# the optimal approximation of the virus in town 2
using homogenous_SIR_model
using Measurements
using DifferentialEquations
#Town Population
N = 10000
I = 1
S = N - I
I_s = 0
R = 0
#Not importing correctly for some reason
function error_squares(model, data)
    err_sum = 0
    for i in 1:(length(data))
        err_sum += (model[i] - data[i])^2
    end

    return err_sum
end

#Infected and severely ill data of the second town
town2_Infected_d27_d80 = [21,29, 25, 30, 28, 34, 28, 54, 57, 92, 73, 80, 109, 102, 128, 135, 163, 150, 211, 196, 233, 247, 283, 286, 332, 371, 390, 404, 467, 529, 598, 641, 704, 702, 788, 856, 854, 955, 995, 1065, 1106, 1159, 1217, 1269, 1298, 1328, 1339, 1383, 1431, 1422, 1414, 1485, 1464, 1480]
town2_Severe_d27_d80 = [3, 3, 4, 7, 3, 8, 7, 5, 9, 13, 15, 3, 20, 13, 11, 20, 16, 11, 15, 18, 27, 24, 28, 36, 41, 35, 41, 55, 63, 66, 72, 80, 90, 104, 109, 115, 127, 135, 147, 162, 163, 186, 194, 200, 216, 223, 241, 249, 258, 275, 277, 299, 302, 300]

#Work back to the assumption that the intervention is implemented at day 37
#Vary beta, epsilon, and p_s

c = 8 #Number of daily contacts on average
gamma = 1/7 # Daily rate of recovery if it takes 7 days to recover typically
gamma_s = 1/14
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month
p_s = measurement(0.2, 0.05)

t_span = (0, 80)
pop0 = [S, I, I_s, R]

Epsilons = range(start = 0.0, step = 0.01, stop = 0.5)
Betas = range(0.02, step = 0.0001, stop = 0.04)
Phis = range(0.2, step = 0.01, stop = 1)


# Store the minimum error as well as the indices it occurs at
# for both infected and severe illnesses
global current_min_error_infected = 10^20
global current_min_indices_infected = (0, 0, 0)
global current_min_error_severe = 10^20
global current_min_indices_severe = (0, 0, 0)


for i in range(1, stop = length(Betas))
    for j in range(1, stop = length(Epsilons))
        for k in range(1, stop = length(Phis))
            local param_int = [c, Betas[i], gamma, alpha, p_s, gamma_s, Epsilons[j], Phis[k]]
            # Create the ODE model of the town
            local model = ODEProblem(town_SIRS_Intervention!, pop0, t_span, param_int)
            local sol = solve(model, saveat=1)
            # The data of interest is the number of infected, obtain from solution as so
            local I_model = [u[2] for u in sol.u]
            local Is_model = [u[3] for u in sol.u]
            local inf_error = error_squares(I_model[27:81], town2_Infected_d27_d80)
            local sev_error = error_squares(Is_model[27:81], town2_Severe_d27_d80)

            # Compare these errors to the current minimum LLSE
            if inf_error < current_min_error_infected
                global current_min_error_infected = inf_error
                global current_min_indices_infected = (i, j, k)
            end
            if sev_error < current_min_error_severe
                global current_min_error_severe = sev_error
                global current_min_indices_severe = (i, j, k)
            end
        

        end
    end
end

println("Beta, epsilon, phi values that minimise infected error: $(Betas[current_min_indices_infected[1]]), $(Epsilons[current_min_indices_infected[2]]), $(Phis[current_min_indices_infected[3]])")
println("Beta, epsilon, phi values that minimise severe illness error: $(Betas[current_min_indices_severe[1]]), $(Epsilons[current_min_indices_severe[2]]), $(Phis[current_min_indices_severe[3]])")
