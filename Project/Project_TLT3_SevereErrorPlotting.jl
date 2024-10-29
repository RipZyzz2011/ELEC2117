using Plots
using DifferentialEquations

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
#Keep beta constant and observe the sensitivity of the severe infection probability
Beta = 0.3513
alpha = 1/30 # Daily rate of resusceptance if the average time for it is a month
