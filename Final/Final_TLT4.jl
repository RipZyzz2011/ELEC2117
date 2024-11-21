using Plots
using DifferentialEquations
using Measurements

#Implements a vaccine at day 30 which makes the population immune for a period of time
function town_SIRSH_Vaccination_decay!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta, gamma, epsilon_s, gamma_s, epsilon_h, gamma_h, alpha, vax, de = param
    #epsilon represents efficacy of intervention
    #phi represents proportion of intervention take-up
    S, I, Is, H, R, V = pop
    if(t < 25)
        vax = 0
        de = 0
    end
    if(t > 60)
        de = 0.2
    end
    #Lambda: Force of Infection
    lambda = c * Beta * (I + Is) 
    dpop[1] = -lambda * S + alpha * R - vax * S #dS
    dpop[2] = lambda * S - gamma * I - epsilon_s * I + lambda*(de) * V#dI
    dpop[3] = epsilon_s * I - gamma_s * Is - epsilon_h * Is #dIs
    dpop[4] = epsilon_h * Is - gamma_h * H #dH
    dpop[5] = gamma * I + gamma_s * Is + gamma_h * H - alpha * R #dR
    dpop[6] = vax * S - lambda*(de) * V #dV
    

end

function town_SIRSH_Vaccination_efficacy!(dpop, pop, param, t)
    N = sum(pop)
    c, Beta, gamma, epsilon_s, gamma_s, epsilon_h, gamma_h, alpha, vax, epsilon_v = param
    #epsilon represents efficacy of intervention
    #phi represents proportion of intervention take-up
    S, I, Is, H, R, V = pop
    if(t < 31)
        vax = 0
        de = 0
    end

    #Lambda: Force of Infection
    lambda = c * Beta * (I + Is) 
    dpop[1] = -lambda * S + alpha * R - vax * S + de * V #dS
    dpop[2] = lambda * S - gamma * I - epsilon_s * I #dI
    dpop[3] = epsilon_s * I - gamma_s * Is - epsilon_h * Is #dIs
    dpop[4] = epsilon_h * Is - gamma_h * H #dH
    dpop[5] = gamma * I + gamma_s * Is + gamma_h * H - alpha * R #dR
    dpop[6] = vax * S - lambda*(de * t/30) * V #dV
    

end

#Q2
#Town Population initial conditions
S = 940
I = 40
Is = 10
H = 5
R = 5
V = 0

N = 1000

#Provided parameters
c = 0.08
beta = 0.0025
gamma = 0.1
epsilon_s = 0.05
gamma_s = 0.05
epsilon_h = 0.1
gamma_h = 0.04
alpha = 0.01
vax = 0.02
de = 0.1

t_span = (0, 200)
pop0 = [S, I, Is, H, R, V]
param = [c, beta, gamma, epsilon_s, gamma_s, epsilon_h, gamma_h, alpha, vax, de]

q3_model = ODEProblem(town_SIRSH_Vaccination_decay!, pop0, t_span, param)
q3_sol = solve(q2_model, saveat = 1)

S_model_q3 = [u[1] for u in q3_sol.u]
I_model_q3 = [u[2] for u in q3_sol.u]
Is_model_q3 = [u[3] for u in q3_sol.u]
H_model_q3 = [u[4] for u in q3_sol.u]
R_model_q3 = [u[5] for u in q3_sol.u]
V_model_q3 = [u[6] for u in q3_sol.u]

N = sum(q2_sol.u[201])
println("$N")

q2_plot = plot()
#plot!(q2_plot, q2_sol.t, S_model_q2, label = "Susceptible")
plot!(q2_plot, q2_sol.t, I_model_q2, label = "Mild Illness")
plot!(q2_plot, q2_sol.t, Is_model_q2, label = "Severe Illness")
plot!(q2_plot, q2_sol.t, H_model_q2, label = "Hospitalized")
#plot!(q2_plot, q2_sol.t, R_model_q2, label = "Recovered")
#plot!(q2_plot, q2_sol.t, V_model_q2, label = "Vaccinated")
plot!(q2_plot, title = "TLT4 SIR Extended Modelling with no Vaccination Decay", xlabel = "Time(Days)", ylabel = "Size of Population Category")

#############################################
#############################################
############################################
#Q3
#Scenario A 
#Town Population initial conditions
S = 940
I = 40
Is = 10
H = 5
R = 5
V = 0

N = 1000

#Provided parameters
c = 0.08
beta = 0.0025
gamma = 0.1
epsilon_s = 0.05
gamma_s = 0.05
epsilon_h = 0.1
gamma_h = 0.04
alpha = 0.01
vax = 0.05
de = 0.05

t_span = (0, 200)
pop0 = [S, I, Is, H, R, V]
param = [c, beta, gamma, epsilon_s, gamma_s, epsilon_h, gamma_h, alpha, vax, de]

q3_model = ODEProblem(town_SIRSH_Vaccination_decay!, pop0, t_span, param)
q3_sol = solve(q3_model, saveat = 1)

S_model_q3 = [u[1] for u in q3_sol.u]
I_model_q3 = [u[2] for u in q3_sol.u]
Is_model_q3 = [u[3] for u in q3_sol.u]
H_model_q3 = [u[4] for u in q3_sol.u]
R_model_q3 = [u[5] for u in q3_sol.u]
V_model_q3 = [u[6] for u in q3_sol.u]

N = sum(q3_sol.u[201])
println("$N")

q3_plot = plot()
plot!(q3_plot, q3_sol.t, S_model_q3, label = "Susceptible")
plot!(q3_plot, q3_sol.t, I_model_q3, label = "Mild Illness")
plot!(q3_plot, q3_sol.t, Is_model_q3, label = "Severe Illness")
plot!(q3_plot, q3_sol.t, H_model_q3, label = "Hospitalized")
plot!(q3_plot, q3_sol.t, R_model_q3, label = "Recovered")
plot!(q3_plot, q3_sol.t, V_model_q3, label = "Vaccinated")
plot!(q3_plot, title = "TLT4 Q3 Scenario A with Vaccination Decay", xlabel = "Time(Days)", ylabel = "Size of Population Category")

#############################################
#############################################
############################################
#
#Q3
#Scenario B
#Town Population initial conditions
S = 940
I = 40
Is = 10
H = 5
R = 5
V = 0

N = 1000

#Provided parameters
c = 0.08
beta = 0.0025
gamma = 0.1
epsilon_s = 0.05
gamma_s = 0.05
epsilon_h = 0.1
gamma_h = 0.04
alpha = 0.01
vax = 0.09
de = 0.1

t_span = (0, 200)
pop0 = [S, I, Is, H, R, V]
param = [c, beta, gamma, epsilon_s, gamma_s, epsilon_h, gamma_h, alpha, vax, de]

q3_model = ODEProblem(town_SIRSH_Vaccination_decay!, pop0, t_span, param)
q3_sol = solve(q3_model, saveat = 1)

S_model_q3 = [u[1] for u in q3_sol.u]
I_model_q3 = [u[2] for u in q3_sol.u]
Is_model_q3 = [u[3] for u in q3_sol.u]
H_model_q3 = [u[4] for u in q3_sol.u]
R_model_q3 = [u[5] for u in q3_sol.u]
V_model_q3 = [u[6] for u in q3_sol.u]

N = sum(q3_sol.u[201])
println("$N")

q3_plot = plot()
plot!(q3_plot, q3_sol.t, S_model_q3, label = "Susceptible")
plot!(q3_plot, q3_sol.t, I_model_q3, label = "Mild Illness")
plot!(q3_plot, q3_sol.t, Is_model_q3, label = "Severe Illness")
plot!(q3_plot, q3_sol.t, H_model_q3, label = "Hospitalized")
plot!(q3_plot, q3_sol.t, R_model_q3, label = "Recovered")
plot!(q3_plot, q3_sol.t, V_model_q3, label = "Vaccinated")
plot!(q3_plot, title = "TLT4 Q3 Scenario B with Vaccination Step-Decay", xlabel = "Time(Days)", ylabel = "Size of Population Category")

