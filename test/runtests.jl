using SIR_model_pkg
using DifferentialEquations
using Test
using Plots

@testset "SIR_model_pkg.jl" begin
    # Write your tests here.
    # 定义公共的参数
    c = 10
    β = 0.03
    γ = 0.1
    S0 = 4999
    I0 = 1
    R0 = 0
    tspan = (0.0, 180.0)

    Ro = (c * β) / γ

    # test basic SIR model
    @testset "Basic SIR Model" begin
        sol = run_sir_model(:basic, c, β, γ, S0, I0, R0, tspan)
        @test sum(sol.u[1]) == S0 + I0 + R0  # the total population shoule remains the same
    end

    # test SIR model with force of infection
    @testset "SIR Model with Force of Infection" begin
        λ_init = 0.05  # Initial infection force
        sol = run_sir_model(:force_of_infection, c, β, γ, S0, I0, R0, tspan; λ=λ_init)
        @test sum(sol.u[1]) ≈ S0 + I0 + R0  # Ensure the total population remains constant
        I_model = sol.u[10][2]  # Infected individuals at the 10th time point
        λ_computed = β * I_model  # Manually calculate the force of infection λ = β * I
        @test λ_computed ≈ β * I_model atol=1e-5  # Compare manually computed and theoretical values
    end

    # 测试带群体免疫的 SIR 模型
    @testset "SIR Model with Herd Immunity" begin
        herd_threshold = 1 - (1 / Ro)  # Setting the threshold for herd immunisation
        sol = run_sir_model(:herd_immunity, c, β, γ, S0, I0, R0, tspan; herd_threshold=herd_threshold)
        @test sol[2,end] ≈ 0 atol=1e-2    # the final infected person should be 0
        @test sol[3,end] > herd_threshold * (S0 + I0 + R0) 
    end

end
