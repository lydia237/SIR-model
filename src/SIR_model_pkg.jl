module SIR_model_pkg

# Import necessary packages
using DifferentialEquations
using Plots
using Revise
using Random

# Export relevant functions for external usage
export sir_model
export run_sir_model
export plot_model_vs_data
export error_num
export optimize_parameters
export plot_overall_model
export plot_error_vs_beta
export run_sensitivity_analysis
export plot_sensitivity_analysis
export plot_infected_and_severe_by_beta  
export analyze_impact_of_coverage

# Include the file containing the model functions
include("model_function.jl")

end
