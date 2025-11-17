#=
    @file beamformer.jl
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
module Beamformer

# Includes, Exports, and Imports
include("./configs.jl")
using .SwellConfigs 

using Match # Module to make nicer looking Enum Match Statements

export process

# Public Function that main.jl calls
function process(method::ProcessingMethod)
    @match method begin 
        $conventional => process_conv()
        $conventional_symlets => process_conv_sym()
        $sbl => process_sbl()
        $sbl_symlets => process_sbl_sym()

        # If not one of the above types, pause
        _ => ErrorException("Please choose a valid Processing type.")
    end
end
    
# Internal Calls
function process_conv()
    println("In the Conventional Processing Function!")
end

function process_conv_sym()
    println("In the Conventional Symlets Processing Function!")
end

function process_sbl()
    println("In the Sparse Bayesian Learning Processing Function!")
end

function process_sbl_sym()
    println("In the Sparse Bayesian Learning with Symlets Processing Function!")
end

end # module