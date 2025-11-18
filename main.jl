#=  @file main.jl
    @author Avik Ghosh
    @date November 16th, 2025

    Purpose: Main location to run each of the numerical routines from. The routines are:
                (1) Conventional Beamforming
                (2) Conventional Beamforming using Symlets
                (3) Sparse Bayesian Learning
                (4) SBL using Symlets

    More detail about each implementation can be found in the implementation specific files
=#

# Includes and Using Statements
include("./src/bf_controller.jl")

using .BFController
using .SwellConfigs

function run()
    # Data is not provided here, but is in a .mat file.
    filePath = joinpath(pwd(), "snapshotData64.mat")

    # Create Config Struct with defaults. Everything else other than the experiment config is default
    # and constant. The experiment config will be updated based on the data file itself.
    T::Type = UInt32
    U::Type = Float64
    V::Type = ComplexF64

    convergence = ConvergenceConfig{T, U}()
    sensor = SensorConfig{U}()
    geometry = GeometryConfig{T, U}()
    sbl = SBLConfig{T}()
    wavelet = WaveletConfig{T}()

    # Data is all in the form of Complex Float64 points, so set experiment config accordingly
    experiment = ExperimentConfig{T, U, V}()

    cfgs = PhysicalConfigs{T, U, V}(convergence, sensor, sbl, wavelet, geometry, experiment)

    # 1. Conventional Beamforming
    method = BFController.conventional
    process(method, cfgs, filePath)

    # # 2. Conventional Beamforming with Symlets
    # method = BFController.conventional_symlets
    # process(method)
    
    # 3. Sparse Bayesian Learning 
    method = BFController.sbl 
    process(method, cfgs, filePath)
    
    # # 4. Sparse Bayesian Learning with Symlets
    # method = BFController.sbl_symlets
    # process(method)
end