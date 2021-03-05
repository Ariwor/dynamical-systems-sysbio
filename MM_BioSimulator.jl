#= A (tutorial) model of a basic enzymatic reaction (E + S <-> C -> E + P) & stochastic simulation (SSA)
 using the BioSimulator.jl package. =#

using BioSimulator, Plots

paths = 1000

# Initialize
model = Network("Michaelis-Menten")

# Define species
model <= Species("S", 301)
model <= Species("E", 130)
model <= Species("C", 0)
model <= Species("P", 0)

# Define reactions
model <= Reaction("dimerization", 0.00166, "S + E --> C")
model <= Reaction("dissociation", 0.0001, "C --> S + E")
model <= Reaction("conversion", 0.1, "C --> P + E")

# Simulate the model
result = simulate(model, Direct(), tfinal = 100.0, save_points = 0:25:100.0)

# Run multiple "simulations" (multiple trajectories)
ensemble = [simulate(model, Direct(), tfinal = 100.0, save_points = 0:25:100.0) for _ in 1:paths]

# Plot (one trajectory)

plot(result[2], summary = :trajectory,
    xlabel = "time", ylabel = "P")

# Plot (mean of multiple trajectories)

plot(ensemble[1], summary = :trajectory,
    xlabel = "time", ylabel = "S",
    label = ["S" "E" "C" "PE"])

plot(ensemble, summary = :mean,
    xlabel = "time", ylabel = "average",
    vars = [1,2], label = ["S" "E"])



    

    
