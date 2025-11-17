#=
    @file configs.jl
    @author Avik Ghosh
    @date November 16th, 2025

    Purpose: Organize all the configs for the Julian implementation of SBL with Wavelet 
                domain processing.
=#

using Parameters
module SBLConfigs

@with_kw struct ConvergenceConfig
    # Convergence Error Thresholds
    error::Float32 = 10e-3
    delay::UInt16 = 200

    # Iteration Parameters
    # These are only used if SBLCfg is configured properly (.iter_flag = True)
    min_iter::UInt16 = 15
    max_iter::UInt32 = 1000
end

@with_kw struct SensorConfig
    # Noise Power Initialization Guess
    noise_power_guess::Float32    = 0.1       # [1/Hz] (Confirm units?)

    # Physical info
    sound_speed::Float32 = 1488                 # [m / s]
end

@with_kw struct GeometryConfig
    # Details about the Array
    num_sensors::UInt16 = 64                    # 64 Sensors is base, can be adjusted as necessary
    sensor_spacing::Float32 = 1.875             # [m]
    geometry::Symbol = :VLA                     # Maybe helpful, but geometry is not something that is used. 
                                                # This is constant/not checked for because our Data is all VLA


    # Known Broken Indices
    broken_indices::Vector{Int32} = Int32[43];

    # Resolution for DoA Estimation
    angles::Vector{Float32} = -90:0.5:90        # [°]
end

@with_kw struct SBLConfig 
    # Number of iterations before you report, ~ Iteration Period
    status_report_num_iters::UInt16 = 150

    # Other Config Parameters
    tic::UInt16 = 1
    iter_flag::Bool = false     # If 1, warns if not converged by max_iter above

    # Fixed Point: Options are {0.5, 2} and 2 is chosen by default
    fixed_point::UInt16 = 2

    # Estimate number of Sources
    num_srcs::UInt16= 3
end

@with_kw struct WaveletConfig
    type::Symbol = :symlet
    order::UInt8 = 8
end

# Struct containing all information for the configs
struct PhysicalConfigs
    Convergence::ConvergenceConfig
    Sensor::SensorConfig
    SBL::SBLConfig
    Wavelet::WaveletConfig
    Geometry::GeometryConfig
end

end # module