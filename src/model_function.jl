# Updated function to define the SIR model with a severe illness compartment (Is) and intervention after 30 days
function sir_model!(dpop, pop, p, t)
    S, I, Is, R = pop  # Susceptible, Infected, Severe illness, Recovered
    c, β, γ, ps, γs, α, epsilon, phi = p  # Add epsilon (efficacy) and phi (coverage) as parameters
    N = S + I + Is + R  # Total population
    
    # Apply intervention after day 30
    if t >= 30.0
        λ = c * (1 - epsilon * phi) * β * I / N  # Adjusted force of infection with intervention
    else
        λ = c * β * I / N  # Original force of infection without intervention
    end

    dpop[1] = -λ * S + α * R  # Change in susceptible
    dpop[2] = λ * S - γ * I  # Change in infected
    dpop[3] = γ * ps * I - γs * Is  # Change in severely ill
    dpop[4] = γ * (1 - ps) * I + γs * Is - α * R  # Change in recovered
end

# Function to run the SIR model with given parameters and initial conditions
function run_sir_model(c, β, γ, ps, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan)
    u0 = [S0, I0, Is0, R0]
    p = [c, β, γ, ps, γs, α, epsilon, phi]  # Include epsilon and phi in parameters
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
function optimize_parameters(beta_range, ps_range, c, γ, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan, infected_days, severe_days, infected_data, severe_data)
    global min_error = Inf  # Initialize the minimum error as infinity
    global best_beta = 0.0  # Initialize the best beta
    global best_ps = 0.0  # Initialize the best ps
    global sol_optimal  # To store the solution with the best parameters

    # Search over all combinations of β and ps
    for β in beta_range
        for ps in ps_range
            sol = run_sir_model(c, β, γ, ps, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan)  # Run the model

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

function plot_error_vs_beta(beta_range, ps, c, γ, γs, α, S0, I0, Is0, R0, tspan, infected_days, severe_days, infected_data, severe_data)
    errors = []
    for β in beta_range
        sol = run_sir_model(c, β, γ, ps, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan)
        infected_model_data = [sol(t)[2] for t in infected_days]
        severe_model_data = [sol(t)[3] for t in severe_days]
        total_error = error_num(infected_model_data, infected_data) + error_num(severe_model_data, severe_data)
        push!(errors, total_error)
    end
    plt = plot(beta_range, errors, xlabel="Beta", ylabel="Error", label="Error as a Function of Beta")
    display(plt)
end

# Function to run sensitivity analysis with varying parameters
function run_sensitivity_analysis(param_sets, S0, I0, Is0, R0, tspan)
    u0 = [S0, I0, Is0, R0]
    solutions = []
    
    for params in param_sets
        prob = ODEProblem(sir_model!, u0, tspan, params)
        sol = solve(prob, saveat=0.1)
        push!(solutions, sol)
    end
    
    return solutions
end

# Function to plot sensitivity analysis with labeled parameter sets and real data points
function plot_sensitivity_analysis(solutions, param_sets, infected_days, infected_data, severe_days, severe_data)
    plt = plot()
    
    # Iterate over solutions and parameters, filter out solutions where the infected population is too large
    for (i, sol) in enumerate(solutions)
        params = param_sets[i]
        β, ps = params[2], params[3]  # Now fetching β and ps for the legend
        
        # Filter out solutions where the final infected population exceeds a threshold (e.g., 100)
        if maximum(sol[2, :]) <= 200  # Set a threshold to filter out extreme solutions
            plot!(plt, sol.t, sol[2, :], label="β=$β, ps=$ps", lw=2)  # Update label with β and ps
        end
    end
    
    # Add real data points as scatter plot on top of the sensitivity analysis
    scatter!(plt, infected_days, infected_data, label="Real Infected Data", color=:red, marker=:circle)

    # Set axis labels
    xlabel!("Time (days)")
    ylabel!("Infected Population (I)")
    
    display(plt)  # Display the plot
end

# Function to plot Infected and Severe Illness populations over time for different β values
function plot_infected_and_severe_by_beta(beta_range, c, γ, ps, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan)
    infected_plot = plot(xlabel="Time (days)", ylabel="Infected Population (I)", legend=:right)
    severe_plot = plot(xlabel="Time (days)", ylabel="Severe Illness Population (Is)", legend=:right)

    # Loop through each beta value and plot the results
    for β in beta_range
        # Run the model for the current beta
        sol = run_sir_model(c, β, γ, ps, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan)
        
        # Plot Infected (I) population over time
        plot!(infected_plot, sol.t, sol[2, :], label="β = $β", lw=2)
        
        # Plot Severe Illness (Is) population over time
        plot!(severe_plot, sol.t, sol[3, :], label="β = $β", lw=2)
    end

    # Display the plots
    display(infected_plot)
    display(severe_plot)
end

# Function to explore impact of coverage and beta on peak infected and severe illness populations
function analyze_impact_of_coverage(beta_range, coverage_range, c, γ, ps, γs, α, S0, I0, Is0, R0, tspan)
    peak_infected = Dict()  # To store peak infected populations
    peak_severe = Dict()    # To store peak severe illness populations

    for β in beta_range
        peak_infected[β] = []
        peak_severe[β] = []

        for φ in coverage_range
            # Run model with current beta and coverage
            sol = run_sir_model(c, β, γ, ps, γs, α, 0.3, φ, S0, I0, Is0, R0, tspan)

            # Record the peak values
            max_infected = maximum(sol[2, :])
            max_severe = maximum(sol[3, :])

            # Append to results for plotting later
            push!(peak_infected[β], max_infected)
            push!(peak_severe[β], max_severe)
        end
    end

    # Plot peak infected population vs coverage for different beta values
    infected_plot = plot(xlabel="Intervention Coverage (ϕ)", ylabel="Peak Infected Population (I)", legend=:topright)
    for β in beta_range
        plot!(infected_plot, coverage_range, peak_infected[β], label="β = $β", lw=2)
    end
    display(infected_plot)

    # Plot peak severe illness population vs coverage for different beta values
    severe_plot = plot(xlabel="Intervention Coverage (ϕ)", ylabel="Peak Severe Illness Population (Is)", legend=:topright)
    for β in beta_range
        plot!(severe_plot, coverage_range, peak_severe[β], label="β = $β", lw=2)
    end
    display(severe_plot)
end

# Initial parameters and values provided by the Department of Health
c = 8  # Average number of contacts per person per day
γ = 0.1429  # Recovery rate (1/7 days)
ps_range = 0.15:0.001:0.25  # Proportion developing severe illness, provided as 15%-25%
γs = 0.0714  # Recovery rate for severe illness (1/14 days)
α = 0.0333  # Re-susceptibility rate (1/30 days)
epsilon = 0.3 # Efficacy of the intervention
phi = 0.8 # Population of the intervention
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

# Define the range of possible values for β and coverage
beta_range = 0.025:0.0001:0.045  # Transmission rate (β) to be optimized
beta_range_modified = 0.025:0.005:0.045
coverage_range = 0.5:0.1:0.9  # Range of intervention coverage values

# Optimize the parameters
best_beta, best_ps, min_error, sol_optimal = optimize_parameters(beta_range, ps_range, c, γ, γs, α, epsilon, phi, S0, I0, Is0, R0, tspan, infected_days, severe_days, infected_data, severe_data)

# Herd immunity threshold
Ro = c * best_beta / γ
pc = 1 - 1/Ro 

# Print the results of the optimization
println("Best β (transmission probability): ", best_beta)
println("Best ps (proportion with severe illness): ", best_ps)
println("Minimum squared error: ", min_error)
println("basic reproduction number: ", Ro)
println("Herd Immunity Threshold: ", pc)

# Define parameter ranges for sensitivity analysis
transmission_probs = [0.025, 0.035, 0.045]  # Varying transmission probability (β)
ps_values = [0.15, 0.2, 0.25]  # Varying proportion with severe illness (ps)

# Combine parameter sets for the simulations
param_sets = [[8, β, ps, 0.1429, 0.0714, 0.0333, epsilon, phi] for β in transmission_probs, ps in ps_values]

# Run sensitivity analysis
solutions = run_sensitivity_analysis(param_sets, S0, I0, Is0, R0, tspan)

# Plot sensitivity analysis with labeled parameter sets and real data points
plot_sensitivity_analysis(solutions, param_sets, infected_days, infected_data, severe_days, severe_data)

# After optimization, run the model over t = 0 to 90 days using the best parameters
full_tspan = (0.0, 210.0)  # Time span for the full simulation
sol_full = run_sir_model(c, best_beta, γ, best_ps, γs, α, epsilon, phi, S0, I0, Is0, R0, full_tspan)

# Plot the overall model with the full time range
plot_overall_model(sol_full, full_tspan)

# Plot the optimal model against the real data
plot_model_vs_data(sol_optimal, infected_days, severe_days, infected_data, severe_data)

# Plot the error against beta
plot_error_vs_beta(beta_range, best_ps, c, γ, γs, α, S0, I0, Is0, R0, tspan, infected_days, severe_days, infected_data, severe_data)

# Run and plot the model for each beta value
plot_infected_and_severe_by_beta(beta_range_modified, c, γ, best_ps, γs, α, epsilon, phi, S0, I0, Is0, R0, full_tspan)

# Run the analysis
analyze_impact_of_coverage(beta_range_modified, coverage_range, c, γ, best_ps, γs, α, S0, I0, Is0, R0, (0.0, 210.0))

