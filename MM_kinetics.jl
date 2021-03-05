#= A model of a basic enzymatic reaction (E + S <-> C -> E + P), numerical simulation
 of the ODEs that describe it and comparison with the Michaelis-Menten approximation. =#

using DifferentialEquations
using Plots; gr()

# ODEs that describe the basic enzymatic reaction E + S <-> C -> E + P
function MMkinetics(dx,x,k,t)
    k_1,k_minus_1,k_2 = k
    dx[1] = k_minus_1*x[3] - k_1*x[1]*x[2] # dS/dt
    dx[2] = (k_minus_1 + k_2)*x[3] - k_1*x[2]*x[1] # dE/dt
    dx[3] = -(k_minus_1 + k_2)*x[3] + k_1*x[2]*x[1] # dC/dt
    dx[4] = k_2*x[3] # dP/dt
end

x0 = [10.0,1.0,0.0,0.0]
k = (1,1,1)
tspan = (0.0,30.0)
MMkinetics_example = ODEProblem(MMkinetics,x0,tspan,k)
sol = solve(MMkinetics_example, Tsit5())

plot(sol,linewidth=1.5,title="Numerical solution of the set of ODEs",
     xaxis="Time (t)",yaxis="Concentration")

# Dynamics with Michaelis-Menten approximation
E_tot = 1
S_m = sol[1,:]
K_m = (k[2] + k[3])/k[1]
C_m = E_tot .* S_m ./(K_m .+ S_m)
E_m = E_tot .- C_m
P_m = x0[1] .- S_m .- C_m

# Graphical comparison with the original numerical simulation
plot(sol.t,[S_m,E_m,C_m,P_m], xscale=:log, xlims = (0.001, 30))
plot!(sol, xscale=:log, xlims = (0.001, 30),line=:dash,title="Numerical solution of the set of ODEs",
     xaxis="Time (t)",yaxis="Concentration")
