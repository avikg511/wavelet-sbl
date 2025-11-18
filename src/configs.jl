#=  @file configs.jl
    @author Avik Ghosh
    @date November 16th, 2025

    Purpose: Organize all the configs for the Julian implementation of SBL with Wavelet 
                domain processing.
=#

using Parameters
module SwellConfigs 

# Exports
export PhysicalConfigs, ConvergenceConfig, SensorConfig, ExperimentConfig, GeometryConfig, SBLConfig, WaveletConfig
export ProcessingMethod
export conventional, conventional_symlets, sbl, sbl_symlets

Base.@kwdef struct ConvergenceConfig{T <: Integer, U <: AbstractFloat}
    # Convergence Error Thresholds
    error::U = 10e-3
    delay::T = 200

    # Iteration Parameters
    # These are only used if SBLCfg is configured properly (.iter_flag = True)
    min_iter::T = 15
    max_iter::T = 1000
end

Base.@kwdef struct SensorConfig{U <: AbstractFloat}
    # Noise Power Initialization Guess
    noise_power_guess::U    = 0.1       # [1/Hz] (Confirm units?)

    # Physical info
    sound_speed::U = 1488.0                 # [m / s]
end

Base.@kwdef struct SBLConfig{T <: Integer}
    # Number of iterations before you report, ~ Iteration Period
    status_report_num_iters::T = 150

    # Other Config Parameters
    tic::T = 1
    iter_flag::Bool = false     # If 1, warns if not converged by max_iter above

    # Fixed Point: Options are {0.5, 2} and 2 is chosen by default
    fixed_point::T = 2

    # Estimate number of Sources
    num_srcs::T = 3
end

Base.@kwdef struct WaveletConfig{T <: Integer}
    type::Symbol = :symlet
    order::T = 8
end

Base.@kwdef struct GeometryConfig{T <: Integer, U <: AbstractFloat}
    # Details about the Array
    num_sensors::T = 64                          # 64 Sensors is base, can be adjusted as necessary
    sensor_spacing::U = 1.875                  # [m]
    geometry::Symbol = :VLA                             # Maybe helpful, but geometry is not something that is used. 
                                                        # This is constant/not checked for because our Data is all VLA

    # Known Broken Indices
    broken_indices::Vector{T} = T[43];

    # Resolution for DoA Estimation
    angles::Vector{U} = -90:0.5:90                      # [°]
end

Base.@kwdef struct ExperimentConfig{T <: Integer, U <: AbstractFloat, V <: Complex}
    # Number of frequencies/tones, sensors, and snapshots/points in time where data is collected
    num_tones::T = 0
    num_sensors::T = 0
    num_snapshots::T = 0

    # For Recording, start/stop time
    start_time::U = 0.0
    stop_time::U = 0.0
    
    # Data, Complex valued Array [num_sensors x num_snapshots x num_tones]
    data::Array{V, 3} = zeros(V, 2, 2, 2)
    all_tones::Vector{U} = U[]
end

# Struct containing all information for the configs
mutable struct PhysicalConfigs{T <: Integer, U <: AbstractFloat, V <: Complex}
    Convergence::ConvergenceConfig{T, U}
    Sensor::SensorConfig{U}
    SBL::SBLConfig{T}
    Wavelet::WaveletConfig{T}
    Geometry::GeometryConfig{T, U}
    Experiment::ExperimentConfig{T, U, V}
end

# Setup for all methods
@enum ProcessingMethod begin
    conventional
    conventional_symlets
    sbl
    sbl_symlets
end

end # module