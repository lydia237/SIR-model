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
    epsilon = 0.3  # Intervention efficacy
    phi = 0.8  # Intervention coverage
    S0 = 5999  # Initial susceptible population
    I0 = 1  # Initial infected population
    Is0 = 0  # Initial severely ill population
    R0 = 0  # Initial recovered population
    tspan = (0.0, 80.0)  # Time span for the simulation

    # Define real data for comparison in plotting tests
    infected_days = 15:80
    severe_days = 21:80
    infected_data = [11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58, 34, 26, 29, 51, 55, 155, 53, 67, 98, 130, 189, 92, 192, 145, 128, 68, 74, 126, 265, 154, 207, 299, 273, 190, 152, 276, 408, 267, 462, 352, 385, 221, 420, 544, 329, 440, 427, 369, 606, 416, 546, 475, 617, 593, 352, 337, 473, 673, 653, 523, 602, 551, 686, 556, 600]  # Infected data for days 15-81
    severe_data = [0, 0, 1, 2, 5, 5, 5, 2, 9, 4, 22, 0, 15, 48, 38, 57, 9, 18, 20, 0, 41, 15, 35, 36, 27, 38, 24, 40, 34, 57, 18, 29, 63, 66, 119, 76, 95, 28, 109, 136, 119, 104, 121, 93, 147, 129, 130, 161, 133, 136, 138, 139, 181, 181, 218, 183, 167, 164, 219, 220]  # Severe illness data for days 21-81


    # Define best β and ps values for optimization testing
    β_optimized = 0.0361  # Optimized transmission rate
    ps_optimized = 0.205  # Optimized proportion with severe illness

    # Test the basic SIR model function
    @testset "Basic SIR Model" begin
        sol = run_sir_model(c, β_optimized, γ, ps_optimized, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan)

        # Ensure the total population remains constant
        total_population = sol[1, :] .+ sol[2, :] .+ sol[3, :] .+ sol[4, :]
        @test all(total_population .≈ S0 + I0 + Is0 + R0)  # Total population should remain unchanged

        # Ensure the model runs correctly without negative values
        @test all(sol[1, :] .>= 0)  # Susceptible should never be negative
        @test all(sol[2, :] .>= 0)  # Infected should never be negative
        @test all(sol[3, :] .>= 0)  # Severe illness should never be negative
        @test all(sol[4, :] .>= 0)  # Recovered should never be negative
    end

    # Test the optimization function
    @testset "Optimization Test" begin
        beta_range = 0.025:0.0001:0.045
        ps_range = 0.15:0.001:0.25
        coverage_range = 0.5:0.1:0.9
        best_beta, best_ps, best_phi, min_error, sol_optimal = optimize_parameters(
            beta_range, ps_range, c, γ, γs, α, epsilon, coverage_range, S0, I0, Is0, R0, tspan, infected_days, severe_days, infected_data, severe_data
        )

        # Check if the optimized parameters are within expected ranges
        @test best_beta in beta_range
        @test best_ps in ps_range
        @test best_phi in coverage_range
        @test min_error >= 0  # Minimum error should be non-negative
    end
end
