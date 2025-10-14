# The value of the concentration of P=HOBr is set to 1 arbitrarily.

#using Pkg
#Pkg.add("Catalyst")
#Pkg.add("OrdinaryDiffEqDefault")
#Pkg.add("Plots")
#Pkg.add("TensorBoardLogger")

using Catalyst
using OrdinaryDiffEqDefault
using Plots
using TensorBoardLogger

function Oregonator_default()
    reaction_network = @reaction_network begin
        @species X(t) Y(t) Z(t)
        @species A P [constant=true] 
        @parameter f
        @parameters k1 k2 k3 k4 k5
        k1, A + Y --> X + P
        k2, X + Y --> 2*P
        k3, A + X --> 2*X + Z
        k4, 2*X --> A + P
        k5, Z --> f*Y
    end
    
    # Simulate the concentration of the different species over time
    variable_concentrations = [:X => 1.0, :Y => 1.0, :Z => 1.0] # initial conditions
    constants = [:k1 => 1.28, :k2 => 2.4*10^6, :k3 => 33.6, :k4 => 2400, :k5 => 1, :A => 0.06, :P => 1, :f => 0.1] # kinetic constants, constant concentrations (values can be adapted for optimal metabolism) and f
    tspan = (0.0, 10.0) # 1500
    
    TensorBoard_logger = TBLogger("logs")
    saved_concentrations = (concentrations, tspan, integrator::ODEProblem)
    saved_values = SavedValues("Float64", "Float64")
    
    log_to_tensorboard = SavingCallback(saved_concentrations, saved_values) # save at every timestep for tensorboard visualization

    oprob = ODEProblem(reaction_network, variable_concentrations, tspan, constants) # ODE problem formulation (can be given a specific algorithm can be determined for optimization), check the constants values
    sol = solve(oprob; callback=log_to_tensorboard, saveat=0.5)
    concentrations = plot(sol)
    savefig(concentrations, "Oregonator_default.pdf")
    return concentrations
end

Oregonator_default()