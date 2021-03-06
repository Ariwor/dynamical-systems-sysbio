#= A  deterministic model of a gene expression network controlled by an ideal AIF controller. 
AIF motif introduced in Briat et al. 2016. Cell Systems. =#

using DifferentialEquations
using Plots; gr()

# ODEs that describe AIF control of a gene expression network

function ideal_AIF_gene_expression(dx,x,k,t)
    mu,kappa,theta,eta,gamma_1,kappa_2,gamma_2 = k
    dx[1] = mu - eta*x[1]*x[2] # dZ1/dt
    dx[2] = theta*x[4] - eta*x[1]*x[2] # dZ2/dt
    dx[3] = kappa*x[1] - gamma_1*x[3] # dX1/dt
    dx[4] = kappa_2*x[3] - gamma_2*x[4] # dX2/dt
end

x0 = [0.0,0.0,0.0,0.0]
k = (3,1,1,50,1,1,2) # Setpoint: mu/theta 
tspan = (0.0,200.0)

ideal_AIF_modeling = ODEProblem(ideal_AIF_gene_expression,x0,tspan,k)
sol = solve(ideal_AIF_modeling, Tsit5())

# Plotting

plot(sol, linewidth=1.5,title="AIF control of a gene expression network",
     xaxis="Time (t)",yaxis="Concentration", label = ["Z1" "Z2" "X1" "X2"])
