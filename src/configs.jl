#=  @file configs.jl
    @author Avik Ghosh
    @date November 16th, 2025

    Purpose: Organize all the configs for the Julian implementation of SBL with Wavelet
                domain processing.
=#

using Parameters
module SwellConfigs

# Export structs and enums
export PhysicalConfigs, ConvergenceConfig, SensorConfig, ExperimentConfig
export GeometryConfig, SBLConfig, WaveletConfig
export ProcessingMethod, conventional, conventional_symlets, sbl, sbl_symlets, sbl_denoising

Base.@kwdef struct ConvergenceConfig{T <: Integer, U <: AbstractFloat}
    # Convergence Error Thresholds
    error::U = 10e-4
    delay::T = 200

    # Iteration Parameters. Bool iter_flag in SBLConfig changes how these are used
    min_iter::T = 15
    max_iter::T = 1000
end

Base.@kwdef struct SensorConfig{U <: AbstractFloat}
    noise_power_guess::U = 0.1      # Initial Guess for Noise Power
    sound_speed::U = 1488.0         # [m/s], in water. Assuming constant with depth
end

Base.@kwdef struct SBLConfig{T <: Integer}
    # Number of iterations before you report, ~ T_{Report Status}
    status_report_num_iters::T = 150

    # Other Config Parameters
    tic::T = 1              # Unused, vestigial from MATLAB Implementation.
    iter_flag::Bool = true  # If True, warns if not converged by max_iter above
    num_srcs::T = 3         # Find this many sources in SBL, throw away other peaks

    # Fixed Point: Options are {0.5, 2} and 2 is chosen by default
    fixed_point::T = 2
end

# Wavelet Configurations - Currently only Symlets supported
Base.@kwdef struct WaveletConfig{T <: Integer}
    type::Symbol = :symlet
    order::T = 8
end

# Geometry Configurations, Defaults based on SwellEx96 VLA
Base.@kwdef struct GeometryConfig{T <: Integer, U <: AbstractFloat}
    # Details about the Array
    num_sensors::T = 64                          # All sensors (with broken ones)
    sensor_spacing::U = 1.875                  # [m] between sensors, vertical
    geometry::Symbol = :VLA                             # Unused, VLA = Vertical Line Array
    angles::Vector{U} = -90:0.5:90                      # [°], Resolution for DoA detection

    # Known Broken Indices - Will be removed while reading MAT data
    broken_indices::Vector{T} = T[43];
end

Base.@kwdef struct ExperimentConfig{T <: Integer, U <: AbstractFloat, V <: Complex}
    # Number of frequencies/tones, sensors, and snapshots/points in time where data is collected
    num_tones::T = 0
    num_sensors::T = 0
    num_snapshots::T = 0

    # For Recording, start/stop time
    start_time::U = 0.0
    stop_time::U = 0.0

    # Data, Complex valued Array [num_sensors x num_snapshots x num_tones]. Dud array now
    data::Array{V, 3} = zeros(V, 2, 2, 2)
    all_tones::Vector{U} = U[]
end

# Struct containing all information for the configs. Passed into almost all functions.
mutable struct PhysicalConfigs{T <: Integer, U <: AbstractFloat, V <: Complex}
    Convergence::ConvergenceConfig{T, U}
    Sensor::SensorConfig{U}
    SBL::SBLConfig{T}
    Wavelet::WaveletConfig{T}
    Geometry::GeometryConfig{T, U}
    Experiment::ExperimentConfig{T, U, V}
end

# Enum to chose processing method at main.jl level
@enum ProcessingMethod begin
    conventional
    conventional_symlets
    sbl
    sbl_symlets
    sbl_denoising
end

end # module
