module SIR_model_pkg

# Write your package code here.
using DifferentialEquations # this pkg is not used in this simple example 
using Plots 
export sir_model 
export run_simulation
export sir_model_herd
export run_sir_model
include("model_function.jl")

end
