using SIR_model_pkg
using DifferentialEquations
using Test
using Plots
using Revise

@testset "SIR_model_pkg.jl" begin
    # Define initial parameters for the tests
    c = 8  # Contact rate
    γ = 0.1429  # Recovery rate (1/7 days)
    ps = 0.2  # Proportion with severe illness
    γs = 0.0714  # Recovery rate for severe illness (1/14 days)
    α = 0.0333  # Re-susceptibility rate (1/30 days)
    S0 = 5999  # Initial susceptible population
    I0 = 1  # Initial infected population
    Is0 = 0  # Initial severely ill population
    R0 = 0  # Initial recovered population
    tspan = (0.0, 25.0)  # Time span for the simulation

    # Define the best β and ps values for optimization testing
    β = 0.0353  # Transmission rate obtained from optimization
    ps_optimized = 0.15  # Proportion with severe illness from optimization

    # Test basic SIR model function
    @testset "Basic SIR Model" begin
        sol = run_sir_model(c, β, γ, ps_optimized, γs, α, S0, I0, Is0, R0, tspan)

        # Ensure the total population remains constant
        total_population = sol[1, :] .+ sol[2, :] .+ sol[3, :] .+ sol[4, :]
        @test all(total_population .≈ S0 + I0 + Is0 + R0)  # Total population should remain unchanged

        # Ensure the model runs correctly without negative values
        @test all(sol[1, :] .>= 0)  # Susceptible should never be negative
        @test all(sol[2, :] .>= 0)  # Infected should never be negative
        @test all(sol[3, :] .>= 0)  # Severe illness should never be negative
        @test all(sol[4, :] .>= 0)  # Recovered should never be negative
    end

end
