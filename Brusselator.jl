# The values of the concentrations of D and E are set to 1 arbitrarily.

using Pkg
Pkg.add("Catalyst")
Pkg.add("OrdinaryDiffEqDefault")
Pkg.add("Plots")
Pkg.add("TensorBoardLogger")

using Catalyst
using OrdinaryDiffEqDefault
using Plots
using TensorBoardLogger

function Brusselator_default()
    reaction_network = @reaction_network begin
        @species X Y
        @parameters A B D E [isconstantspecies=true] 
        @parameters k1 k2 k3 k4
        
        k1, A --> X
        k2, B + X --> Y + D
        k3, 2*X + Y --> 3*X
        k4, X --> E
    end

    # Simulate the concentration of the different species over time
    variable_concentrations = [:X => 1.0, :Y => 1.0] # initial conditions
    constants = [:k1 => 1.28, :k2 => 2.4*10^6, :k3 => 33.6, :k4 => 2400.0, :A => 0.06, :B => 0.02, :D => 1.0, :E => 1.0] # kinetic constants, constant concentrations (values to be adapted for optimal metabolism)
    tspan = (0.0, 10.0) # 100
    
    oprob = ODEProblem(reaction_network, variable_concentrations, tspan, constants) # ODE problem formulation (can be given a specific algorithm can be determined for optimization), check the constants values
    
    #save at every timestep for tensorboard visualization

    sol = solve(oprob)
    concentrations = plot(sol)
    savefig(concentrations, "Brusselator_default.pdf")
    return concentrations
end

Brusselator_default()