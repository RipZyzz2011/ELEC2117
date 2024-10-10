
using Pkg
Pkg.add(path="C:/Users/hamis/.julia/dev/homogenous_SIR_model/")
using homogenous_SIR_model

# Run the model with some initial conditions
pop0 = [4999, 1, 0]
tspan = (0.0, 50.0)
c = 10
Beta_c = 0.03
gamma = 0.1
p_c = 0
param = [c, Beta_c, gamma, p_c] #c, Beta_c, gamma

model = define_town_model(:herd, param, pop0, tspan)
solved = solve_system(model)
plot_model_solution(solved)