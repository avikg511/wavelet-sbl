#=
    @file main.jl
    @author Avik Ghosh
    @date November 16th, 2025

    Purpose: Main location to run each of the numerical routines from. The routines are:
                (1) Conventional Beamforming
                (2) Conventional Beamforming using Symlets
                (3) Sparse Bayesian Learning
                (4) SBL using Symlets

            More detail about each implementation can be found in the implementation specific files
=#
module MainModule

# Includes and Using Statements
include("./src/beamformer.jl")

using .Beamformer

# 1. Conventional Beamforming
method = Beamformer.conventional
process(method)

# # 2. Conventional Beamforming with Symlets
# method = Beamformer.conventional_symlets
# process(method)
# 
# # 3. Sparse Bayesian Learning 
# method = Beamformer.sbl 
# process(method)
# 
# # 4. Sparse Bayesian Learning with Symlets
# method = Beamformer.sbl_symlets
# process(method)

end # module