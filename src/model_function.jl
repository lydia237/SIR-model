# basic model
function sir_model!(dpop, pop, p, t)
    S, I, R = pop
    c, β, γ = p
    N = S + I + R
    λ = c * β * I / N

    dpop[1] = -λ * S  # Susceptible
    dpop[2] = λ * S - γ * I # Infected
    dpop[3] = γ * I # Recovered
end

# force of infection
function sir_model!(dpop, pop, p, t, λ)
    S, I, R = pop
    c, β, γ = p

    dpop[1] = -λ * S
    dpop[2] = λ * S - γ * I
    dpop[3] = γ * I
end

# herd immunity threshold
function sir_model_herd!(dpop, pop, p, t, herd_threshold)
    S, I, R = pop
    c, β, γ = p
    N = S + I + R

    # If the recovered person exceeds the herd immunity threshold
    # the transmission rate β becomes 0
    if R >= herd_threshold * N
        β = 0
    end

    λ = c * β * I / N

    dpop[1] = -λ * S
    dpop[2] = λ * S - γ * I
    dpop[3] = γ * I
end

# drive function
function run_sir_model(model_type::Symbol, c, β, γ, S0, I0, R0, tspan; λ=nothing, herd_threshold=nothing)
    u0 = [S0, I0, R0]
    p = [c, β, γ]

    if model_type == :basic
        prob = ODEProblem(sir_model!, u0, tspan, p)
    elseif model_type == :force_of_infection && λ !== nothing
        prob = ODEProblem((dpop, pop, p, t) -> sir_model!(dpop, pop, p, t, λ), u0, tspan, p)
    elseif model_type == :herd_immunity && herd_threshold !== nothing
        prob = ODEProblem((dpop, pop, p, t) -> sir_model_herd!(dpop, pop, p, t, herd_threshold), u0, tspan, p)
    #=
    else
        error("Unknown or incomplete model type")
    =#
    end
    sol = solve(prob)  # Assign the solution to `sol`
    
    # Plotting the result
    return plot(sol, title="SIR Model", label=["Susceptible" "Infected" "Recovered"], lw=2, legend=:bottomright)
end

# basic run function 
function run_simulation(c, β, γ, S0, I0, R0, tspan)
    u0 = [S0, I0, R0]
    p = [c, β, γ]
    prob = ODEProblem(sir_model!, u0, tspan, p)
    sol = solve(prob, saveat = 1.0) 
    I_model = sol.u[31][2] #infected population at t = 30
    R_model = sol.u[31][3] #recoverd population at t = 30
    I_proportion = (I_model./10000).*100
    R_proportion = (R_model./10000).*100
    return I_proportion + R_proportion
end


# Function to calculate peak infection
function calculate_peak_infection(sol)
    infected_population = sol[2, :]
    peak_infected = maximum(infected_population)
    peak_time = sol.t[argmax(infected_population)]
    return peak_infected, peak_time
end

function run_simulation_peak(c, β, γ, S0, I0, R0, tspan)
    u0 = [S0, I0, R0]
    p = [c, β, γ]
    prob = ODEProblem(sir_model!, u0, tspan, p)
    sol = solve(prob, saveat = 1)
    return sol
end

4.2
function calculate_herd_immunity_threshold(c, β, γ)
    Ro = c * β / γ
    herd_immunity_threshold = 1 - (1 / Ro)
    return herd_immunity_threshold
end
#=4.2
# Example using the refined β value
c = 15
β_refined = 0.046  # Example refined β from model fitting
γ = 0.1  # Recovery rate (as given in your model)

# Calculate herd immunity threshold
herd_threshold = calculate_herd_immunity_threshold(c, β_refined, γ)
println("Herd Immunity Threshold: $herd_threshold")
=#
#=
4.1
days = [1:1:30]
infected_data = [5, 10, 19, 37, 71, 136, 260, 486, 882, 1516, 2399, 3407, 4300, 4882, 5116, 5080, 4875, 4582, 4251, 3913, 3583, 3271,
2979, 2708, 2460, 2233, 2026, 1837, 1665, 1509]
sol = run_simulation_peak(15, 0.048, 0.1, 9995, 5, 0, (0.0, 30.0))
function error_num(model_data, data)
    if length(model_data) != length(data)
        println("error")
    end
    return sum((model_data[i] - data[i])^2 for i in eachindex(model_data, data))
end
infected_model_data = [sol(t)[2] for t in 1:1:30]
squared_error = error_num(infected_model_data, infected_data)
println("Squared error: ", squared_error)
plot(sol.t, sol[2, :], label="Modeled Infected (u[2])", xlabel="Time (days)", ylabel="Number Infected", title="Outbreak vs Modeled Data")
scatter!(days, infected_data, label="Real Data", color=:red, marker=:circle)
=#

#=3.2
#Example of calculating peak infection
sol = run_simulation_peak(15, 0.03, 0.1, 9995, 5, 0, (0.0, 30.0))
peak_infected, peak_time = calculate_peak_infection(sol)
println("Peak Infection: $peak_infected at time $peak_time")
=#

#2 run_sir_model(:basic, 15, 0.03, 0.1, 9995, 5, 0, (0.0, 30.0))
#3.1 run_simulation(15, 0.03, 0.1, 9995, 5, 0, (0.0, 30.0))