module SIR_model_pkg

# Import necessary packages
using DifferentialEquations
using Plots

# Export relevant functions for external usage
export sir_model
export run_sir_model
export plot_model_vs_data
export error_num
export optimize_parameters
export plot_overall_model

# Include the file containing the model functions
include("model_function.jl")

end
