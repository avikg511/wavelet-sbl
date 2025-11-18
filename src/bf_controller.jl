#=  @file beamformer.jl
    @author Avik Ghosh
    @date November 16th, 2025

    Purpose: main.jl should not deal with a lot of configuration data. The user should just be 
            able to simply configure the code from an upper level and have it work. This file
            will process the configuration, dispatching, etc. based on configs from main.jl

        NOTE: Julia's module import system is a bit finnicky. If one module is imported in two 
            different locations, they are treated as different types. This means each module must
            be imported at one location. This means beamformer takes care of the actual imports and
            main.jl actually is a 'child' of beamformer.jl. This file is the upper level
=# 

# Includes
include("./configs.jl")

module BFController

# Imports
using ..SwellConfigs 
using Match             # Module to make nicer looking Enum Match Statements
using MAT               # For reading .mat files

# Individual Beamforming Modules Per Algorithm
include("./conventional.jl")
using .ConventionalBeamforming

# Export public process function
export process

# Public Function that main.jl calls
function process(method::ProcessingMethod, all_cfgs::PhysicalConfigs{T, U, V}, filePath::String) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    # First, overwrite the Experiment Config with data from the file
    all_cfgs.Experiment = process_data(filePath, all_cfgs.Geometry, T, U, V)

    # Now, based on the method chosen, dispatch to the correct processing routine
    @match method begin 
        $conventional => process_conv(all_cfgs)
        $conventional_symlets => process_conv_sym(all_cfgs)
        $sbl => process_sbl(all_cfgs)
        $sbl_symlets => process_sbl_sym(all_cfgs)

        # If not one of the above types, throw an error
        _ => throw(ErrorException("Please choose a valid Processing type."))
    end
end

# Purpose: Reads an input data file and returns necessary data/measurements for
#     beamforming. This is based on a known structure to the SwellEx96 Dataset
function readDataFile(filePath, brokenSensors, ::Type{T}, ::Type{U}, ::Type{V}) where {T <: Integer, U <:AbstractFloat, V <: Complex}
   # Import .mat file processing. Contains multiple variables
   dataFile = matread(filePath);

   # Extract all the fields/values from the dataFile
   dataMat = dataFile["snapshotsAmplitude"];

   # Delete sensor rows for broken sensors
   for i in brokenSensors 
      dataMat = dataMat[1:end .!= i, :, :]
   end

   # Get all the information about the data and return. Get all_tones first for num_tones
   all_tones = vec(dataFile["allTones"])

   # Return all the information about the data in the form of an ExperimentConfig Struct
   return ExperimentConfig{T, U, V}(
        num_tones = convert(T, size(all_tones)[1]),
        num_sensors = convert(T, size(dataMat)[1]),
        num_snapshots = convert(T, size(dataMat)[2]),

        start_time = convert(U, dataFile["startTime"]), 
        stop_time = convert(U, dataFile["stopTime"]),
        data = convert(Array{V, 3}, dataMat),
        all_tones = convert(Vector{U}, all_tones)
   )
end

# This function does all the templating for the data extraction, though defaults to 
#   UInt32, Float64, ComplexF64 for reasonable precision
function process_data(mat_path::String, geometry::GeometryConfig, ::Type{T}=UInt32, ::Type{U}=Float64, ::Type{V}=ComplexF64) where {T <: Integer, U <:AbstractFloat, V <: Complex}
    # Read data from file, removing broken sensors
    exp_info = readDataFile(mat_path, geometry.broken_indices, T, U, V)
    return exp_info
end
    
# Internal Calls
function process_conv(all_cfgs::PhysicalConfigs{T, U, V}) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    ConventionalBeamforming.conventional_bf(all_cfgs)
end

function process_conv_sym(all_cfgs::PhysicalConfigs)
    println("In the Conventional Symlets Processing Function!")
end

function process_sbl(all_cfgs::PhysicalConfigs)
    println("In the Sparse Bayesian Learning Processing Function!")
end

function process_sbl_sym(all_cfgs::PhysicalConfigs)
    println("In the Sparse Bayesian Learning with Symlets Processing Function!")
end

end # module