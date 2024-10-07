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
    else
        error("Unknown or incomplete model type")
    end

    return solve(prob)
end

# basic run function 
function run_simulation(c, β, γ, S0, I0, R0, tspan)
    u0 = [S0, I0, R0]
    p = [c, β, γ]
    prob = ODEProblem(sir_model, u0, tspan, p)
    sol = solve(prob, saveat = 0.2)
    return sol
end