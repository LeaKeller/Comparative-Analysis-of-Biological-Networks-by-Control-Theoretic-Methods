# The value of the concentration of P=HOBr is set to 1 arbitrarily.

#using Pkg
#Pkg.add("Catalyst")
#Pkg.add("OrdinaryDiffEqDefault")
#Pkg.add("Plots")
#Pkg.add("TensorBoardLogger")

using Catalyst
using OrdinaryDiffEqDefault
@time using Plots
using TensorBoardLogger

function Gyorgyri_default()
    reaction_network = @reaction_network begin
        @species X(t) Y(t) Z(t) V(t)
        @species A M H C [constant=true] 
        @parameters k1 k2 k3 k4 k5 k6 k7
        k1, Y + X + H --> 2*V
        k2, Y + A + 2*H --> V + X
        k3, 2*X --> V
        k4, (1/2)*X + A + H --> X + Z
        k5, X + Z --> (1/2)*X
        k6, V+Z --> Y
        k7, Z + M --> 0
    end
    
    # Simulate the concentration of the different species over time
    variable_concentrations = [:X => 1.0, :Y => 1.0, :Z => 1.0, :V => 1.0] 
    constants = [:M => 0.25, :A => 0.1, :H => 0.26, :C => 0.000833, :k1 => 4.0e6*0.25^(-2), :k2 => 2.0*0.25^(-3), :k3 => 3000*0.25^(-1), :k4 => 55.2*0.25^(-5/2), :k5 => 7000*0.25^(-1), :k6 => 0.09*0.25^(-1), :k7 => 0.23*0.25^(-1)] # constant concentrations (values can be adapted for optimal metabolism) and kinetic contants
    tspan = (0.0, 10.0) # 100
    
    oprob = ODEProblem(reaction_network, variable_concentrations, tspan, constants; combinatoric_ratelaws=false) # ODE problem formulation (can be given a specific algorithm can be determined for optimization), check the constants values
    
    #save at every timestep for tensorboard visualization

    sol = solve(oprob)
    concentrations = @time plot(sol)
    savefig(concentrations, "Gyorgyri_default.pdf")
    return concentrations
end

function Gyorgyri_stochastic()
    reaction_network = @reaction_network begin
        @species X Y Z V
        @parameters A M H C [isconstantspecies=true] 
        @parameters k1 k2 k3 k4 k5 k6 k7
        k1, Y + X + H --> 2*V
        k2, Y + A + 2*H --> V + X
        k3, 2*X --> V
        k4, (1/2)*X + A + H --> X + Z
        k5, X + Z --> (1/2)*X
        k6, V+Z --> Y
        k7, Z + M --> 0
    end
    
    # Simulate the concentration of the different species over time
    variable_concentrations = [[:X => 1.0, :Y => 1.0, :Z => 1.0, :V => 1.0], [:X => 0.4, :Y => 1.2, :Z => 0.5, :V => 0.5], [:X => 0.8, :Y => 1.1, :Z => 0.3, :V => 0.5]] # variation of initial conditions (trials that probaly makes no sense biochemically)
    constants = [:M => 0.25, :A => 0.1, :H => 0.26, :C => 0.000833, :k1 => 4.0e6*0.25^(-2), :k2 => 2.0*0.25^(-3), :k3 => 3000*0.25^(-1), :k4 => 55.2*0.25^(-5/2), :k5 => 7000*0.25^(-1), :k6 => 0.09*0.25^(-1), :k7 => 0.23*0.25^(-1)] # constant concentrations (values can be adapted for optimal metabolism) and kinetic contants
    tspan = (0.0, 10.0) # 100

    # define the concentrations plots : concentrations = ...
    for i = 1:3
        oprob = ODEProblem(reaction_network, variable_concentrations[i], tspan, constants; combinatoric_ratelaws=false) # ODE problem formulation (can be given a specific algorithm can be determined for optimization), check the constants values
        println(oprob)
        sol = solve(oprob)
        concentrations = @time plot(sol)
    end

    #save at every timestep for tensorboard visualization ...

    savefig(concentrations, "Gyorgyri_stochastic.pdf")
    return concentrations
end

Gyorgyri_default()
Gyorgyri_stochastic()