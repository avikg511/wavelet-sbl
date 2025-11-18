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

function run()
    # Data is not provided here, but is in a .mat file.
    filePath = joinpath(pwd(), "snapshotData64.mat")

    # Create Config Struct with defaults. Everything else other than the experiment config is default
    # and constant. The experiment config will be updated based on the data file itself.
    T::Type = UInt32
    U::Type = Float64
    V::Type = Complex{U}            # Complex type based on U, used for complex pressures in data

    convergence     = BFController.SwellConfigs.ConvergenceConfig{T, U}()
    sensor          = BFController.SwellConfigs.SensorConfig{U}()
    geometry        = BFController.SwellConfigs.GeometryConfig{T, U}()
    sbl             = BFController.SwellConfigs.SBLConfig{T}()
    wavelet         = BFController.SwellConfigs.WaveletConfig{T}()
    experiment      = BFController.SwellConfigs.ExperimentConfig{T, U, V}()

    cfgs = BFController.SwellConfigs.PhysicalConfigs{T, U, V}(convergence, sensor, sbl, wavelet, geometry, experiment)

    # 1. Conventional Beamforming
    method = BFController.conventional
    BFController.process(method, cfgs, filePath)

    # # 2. Conventional Beamforming with Symlets
    # method = BFController.conventional_symlets
    # BFController.process(method, cfgs, filePath)

    # 3. Sparse Bayesian Learning 
    # method = BFController.sbl 
    # BFController.process(method, cfgs, filePath)
    
    # # 4. Sparse Bayesian Learning with Symlets
    # method = BFController.sbl_symlets
    # BFController.process(method, cfgs, filePath)
end