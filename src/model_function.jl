# Function to define the SIR model with a severe illness compartment (Is)
function sir_model!(dpop, pop, p, t)
    S, I, Is, R = pop  # Susceptible, Infected, Severe illness, Recovered
    c, β, γ, ps , γs, α= p  # Parameters: contact rate, transmission rate, recovery rate, proportion with severe illness, recovery rate for severe, re-susceptibility rate
    N = S + I + Is + R  # Total population
    λ = c * β * I / N  # Force of infection

    dpop[1] = -λ * S + α * R  # Change in susceptible
    dpop[2] = λ * S - γ * I  # Change in infected
    dpop[3] = γ * ps * I - γs * Is  # Change in severely ill
    dpop[4] = γ * (1 - ps) * I + γs * Is - α * R  # Change in recovered
end

# Function to run the SIR model with given parameters and initial conditions
function run_sir_model(c, β, γ, ps , γs, α, S0, I0, Is0, R0, tspan)
    u0 = [S0, I0, Is0, R0]
    p = [c, β, γ, ps , γs, α]
    prob = ODEProblem(sir_model!, u0, tspan, p)
    sol = solve(prob, saveat=0.1)  # Save data at 0.1 time steps for smoother plotting
    return sol
end

# Function to calculate and plot the model predictions and compare with real data
function plot_model_vs_data(sol, infected_days, severe_days, infected_data, severe_data)
    plt = plot()
    
    # Plot model output for infected and severe illness
    plot!(plt, sol.t, sol[2, :], label="Modeled Infected", xlabel="Time (days)", ylabel="Population", lw=2)
    plot!(plt, sol.t, sol[3, :], label="Modeled Severe Illness", lw=2)
    
    # Scatter plot for real data
    scatter!(plt, infected_days, infected_data, label="Real Infected Data", color=:red, marker=:circle)
    scatter!(plt, severe_days, severe_data, label="Real Severe Illness Data", color=:blue, marker=:square)
    
    display(plt)  # Display the plot
end

# Function to calculate the error between model predictions and real data
function error_num(model_data, data)
    return sum((model_data[i] - data[i])^2 for i in eachindex(model_data, data))  # Sum of squared errors
end

# Function to optimize β (transmission probability) and ps (proportion of severe illness)
function optimize_parameters(beta_range, ps_range, c, γ, γs, α, S0, I0, Is0, R0, tspan, infected_days, severe_days, infected_data, severe_data)
    global min_error = Inf  # Initialize the minimum error as infinity
    global best_beta = 0.0  # Initialize the best beta
    global best_ps = 0.0  # Initialize the best ps
    global sol_optimal  # To store the solution with the best parameters

    # Search over all combinations of β and ps
    for β in beta_range
        for ps in ps_range
            sol = run_sir_model(c, β, γ, ps, γs, α, S0, I0, Is0, R0, tspan)  # Run the model

            # Extract model predictions for infected and severe illness on the specified days
            infected_model_data = [sol(t)[2] for t in infected_days]
            severe_model_data = [sol(t)[3] for t in severe_days]

            # Calculate the total error (sum of infected and severe illness errors)
            total_error = error_num(infected_model_data, infected_data) + error_num(severe_model_data, severe_data)

            # Update the best parameters if this combination gives a smaller error
            if total_error < min_error
                global min_error = total_error
                global best_beta = β
                global best_ps = ps
                global sol_optimal = sol  # Save the best solution
            end
        end
    end
    
    return best_beta, best_ps, min_error, sol_optimal  # Return the best parameters and the optimal solution
end

function plot_overall_model(sol, tspan)
    plt = plot()
    
    # Plot the four compartments: Susceptible, Infected, Severe illness, Recovered
    plot!(plt, sol.t, sol[1, :], label="Susceptible (S)", xlabel="Time (days)", ylabel="Population", lw=2, legend=:topright)
    plot!(plt, sol.t, sol[2, :], label="Infected (I)", lw=2)
    plot!(plt, sol.t, sol[3, :], label="Severe Illness (Is)", lw=2)
    plot!(plt, sol.t, sol[4, :], label="Recovered (R)", lw=2)
    
    display(plt)  # Display the plot
end

# Initial parameters and values provided by the Department of Health
c = 8  # Average number of contacts per person per day
γ = 0.1429  # Recovery rate (1/7 days)
ps_range = 0.15:0.001:0.25  # Proportion developing severe illness, provided as 15%-25%
γs = 0.0714  # Recovery rate for severe illness (1/14 days)
α = 0.0333  # Re-susceptibility rate (1/30 days)
S0 = 5999  # Initial susceptible population
I0 = 1  # Initial infected population
Is0 = 0  # Initial severely ill population
R0 = 0  # Initial recovered population
tspan = (0.0, 25.0)  # Time span in days

# Real data provided for optimization
infected_days = 15:25
severe_days = 21:25
infected_data = [11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58]  # Infected data for days 15-25
severe_data = [0, 0, 1, 2, 5]  # Severe illness data for days 21-25

# Define the range of possible values for β
beta_range = 0.010:0.0001:0.050  # Transmission rate (β) to be optimized

# Optimize the parameters
best_beta, best_ps, min_error, sol_optimal = optimize_parameters(beta_range, ps_range, c, γ, γs, α, S0, I0, Is0, R0, tspan, infected_days, severe_days, infected_data, severe_data)

# Print the results of the optimization
println("Best β (transmission probability): ", best_beta)
println("Best ps (proportion with severe illness): ", best_ps)
println("Minimum squared error: ", min_error)

# After optimization, run the model over t = 0 to 90 days using the best parameters
full_tspan = (0.0, 210.0)  # Time span for the full simulation
sol_full = run_sir_model(c, best_beta, γ, best_ps, γs, α, S0, I0, Is0, R0, full_tspan)

# Plot the overall model with the full time range
plot_overall_model(sol_full, full_tspan)

# Plot the optimal model against the real data
plot_model_vs_data(sol_optimal, infected_days, severe_days, infected_data, severe_data)