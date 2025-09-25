# The value of the concentration of P=HOBr is set to 1 arbitrarily.

using Pkg
Pkg.add("Catalyst")
Pkg.add("OrdinaryDiffEqDefault")
Pkg.add("Plots")
Pkg.add("TensorBoardLogger")

using Catalyst
using OrdinaryDiffEqDefault
using Plots
using TensorBoardLogger

function Gyorgyri_default()
    reaction_network = @reaction_network begin
        k1, Y + X + H --> 2*V
        k2, Y + A + 2*H --> V + X
        k3, 2*X --> V
        k4, (1/2)*X + A + H --> X + Z
        k5, X + Z --> (1/2)*X
        k6, V+Z --> Y
        k7, Z + M --> 0
    end
    
    # Simulate the concentration of the different species over time
    variable_concentrations = [:X => 1.0, :Y => 1.0, :Z => 1.0, :V => 1.0] # initial conditions
    constants = [:M => 0.25, :A => 0.1, :H => 0.26, :C => 0.000833, :k1 => 4.0e6*0.25^(-2), :k2 => 2.0*0.25^(-3), :k3 => 3000*0.25^(-1), :k4 => 55.2*0.25^(-5/2), :k5 => 7000*0.25^(-1), :k6 => 0.09*0.25^(-1), :k7 => 0.23*0.25^(-1), :A => 0.06, :P => 1, :f => 0.1] # constant concentrations (values can be adapted for optimal metabolism) and kinetic contants
    tspan = (0.0, 10.0) # 100
    
    oprob = ODEProblem(reaction_network, variable_concentrations, tspan, constants; combinatoric_ratelaws=false) # ODE problem formulation (can be given a specific algorithm can be determined for optimization), check the constants values
    
    ... #save at every timestep for tensorboard visualization
    
    sol = solve(oprob)
    concentrations = plot(sol)
    savefig(concentrations, "Gyorgyri_default.pdf")
    return concentrations
end

Gyorgyri_default()