#=  @file sbl.jl
    @author Avik Ghosh
    @date November 18th, 2025
    Purpose: Implements the Sparse Bayesian Learning algorithm for 
                beamforming applications using the configs from configs.jl
=#

module SparseBayesianLearning

# Imports
using ..SwellConfigs
using Plots
using LinearAlgebra
using Peaks

function sbl(all_cfgs::PhysicalConfigs{T, U, V}) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    # Create matrices
    num_tones = all_cfgs.Experiment.num_tones
    num_sensors = all_cfgs.Geometry.num_sensors - length(all_cfgs.Geometry.broken_indices)
    num_angles = length(all_cfgs.Geometry.angles)

    # A is the frequency augmented sensing matrix
    A = zeros(Complex{U}, num_tones, num_sensors, num_angles)

    # Generate Sensing Matrices for each frequency and embed into matrix A
    for tone_i in 1:num_tones
        generate_sensing_matrix!(@view(A[tone_i, :, :]), all_cfgs, tone_i)
    end

    report = sparse_bayesian_learning(A, all_cfgs);

    # Sample Plot of final γ values
    plt = plot(report["final_iter_γ"], xlabel="Angle Index", ylabel="Final γ Values", title="SBL Final γ Values per Angle")
    display(plt)
    savefig(pwd() * "/angle_vs_gamma.png")
    return report
end

function generate_sensing_matrix!(sensing_mat::SubArray{Complex{U}, 2}, all_cfgs::PhysicalConfigs{T, U, V}, tone_index::Int) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    # Extract wavelength, which is just speed divided by the tone transmited, λ = c / f
    λ_transmit = all_cfgs.Sensor.sound_speed / all_cfgs.Experiment.all_tones[tone_index];
    
    # Sensor spacing divided by the wavelength itself will be the phase shift of an incoming plane wave between adjacent sensors
    ϕ_shift = 2 * pi * all_cfgs.Geometry.sensor_spacing / λ_transmit;
    N = all_cfgs.Geometry.num_sensors - length(all_cfgs.Geometry.broken_indices);
    sin_θs = sind.(all_cfgs.Geometry.angles);

    # Make a vector from 1 : N but remove elements in all_cfgs.Geometry.broken_indices
    # Then subtract 1 to make it zero indexed. TODO: This is potentially buggy in case
    #  broken indices is zero indexed already. Just following MATLAB implementation here.
    sensor_inds = setdiff(1:all_cfgs.Geometry.num_sensors, all_cfgs.Geometry.broken_indices) .- 1

    # Sensing Matrix
    sensing_mat .= 1 / √(N) * exp.( 1im * ϕ_shift * (reshape(sensor_inds, :, 1) * reshape(sin_θs', 1, :)));
end

function sparse_bayesian_learning(A::Array{Complex{U}, 3}, all_cfgs::PhysicalConfigs{T, U, V}) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    println("In the SBL function!")
    
    # Confirm data is properly set up
    @assert ndims(all_cfgs.Experiment.data) == 3 "Data must be a 3D Array of [sensors, snapshots, freqs]"
    
    # Reorganize Data from [sensors, snapshots, freqs] to [freqs, sensors, snapshots]
    data = permutedims(all_cfgs.Experiment.data, [3, 1, 2])

    # Extract all the configs again
    num_tones = all_cfgs.Experiment.num_tones
    num_sensors = all_cfgs.Geometry.num_sensors - length(all_cfgs.Geometry.broken_indices)
    num_angles = length(all_cfgs.Geometry.angles)
    num_sources = all_cfgs.SBL.num_srcs
    num_snapshots = all_cfgs.Experiment.num_snapshots

    # Initializations - noise power, posterior on x, and l1 error vector
    σ_c = ones(num_tones) * all_cfgs.Sensor.noise_power_guess    
    x_post = zeros(Complex{U}, num_tones, num_angles, num_snapshots)
    l1_errors = zeros(all_cfgs.Convergence.max_iter)

    # Convergence/Error Iterations
    Δγ_max = Inf
    num_iters = 0                   # Only updated upon convergence or leaving the loop
    iteration_L1 = 0                # Iteration where we get min loss
    γmin = zeros(U, num_angles)             # Minimum gamma over iterations
  
    # Set up γ vector and numerator/denominator matrices
    γ = ones(U, num_angles)
    γ_numer = zeros(U, num_tones, num_angles)
    γ_denom = zeros(U, num_tones, num_angles)

    # Sample Covariance Matrix Σ over sensors per frequency/tone
    Σ_sensors = zeros(V, num_tones, num_sensors, num_sensors)
     
    # Get covariance matrices for each tone
    for tone_i in 1:num_tones
        # Get data matrix at each frequency, and set up the base covariance per sensor
        data_i = data[tone_i, :, :]
        Σ_sensors[tone_i, :, :] = (data_i * data_i') / num_snapshots
    end

    ## START OF ITERATIONS ## 
    for iter_i in 1:all_cfgs.Convergence.max_iter
        γ_prev = copy(γ)

        # Gamma update loop, does it in place for each frequency
        gamma_update!(data, A, γ_numer, γ_denom, σ_c, vec(γ), num_tones, num_sensors, num_snapshots)

        # Fixed point equation for gamma
        γ = γ .* transpose( (sum(γ_numer, dims=1) ./ (sum(γ_denom, dims=1))).^(1 / all_cfgs.SBL.fixed_point))
        
        # σ and L2 error using unbiased posterior update
        # Unlike MATLAB, findpeaks does not take in the SORTSTR or NPEAKS commands. We need to 
        #   sort manually.
        pks = findpeaks(vec(γ))

        locs = pks.indices
        vals = pks.heights

        ## TODO: Peaks.jl allows for prominence as well (so separating between nearby peaks, etc. 
        #       so this should be implemented later on)

        # Now we just get the values to sort by (max peak down/descending) and the indices. Get the sorting 
        #   permutation from the sources themselves and keep only the first num_sources of them
        # @show length( loddcs ) length( vals ) 
        @assert length(locs) > num_sources "There are less sources found than peaks from the findmaxima function. Fix." 
        val_permute = sortperm(vals, rev = true)[1:num_sources]
        source_locs = locs[val_permute]

        # Extract A for the peaks, basically index over the angles to find what peaks are reasonable
        #   given num_sources, and apply that to the actual A matrix and extract data for those angles only
        A_peak = A[:, :, source_locs]

        for tone_i in 1:num_tones
            # Get the active replicas for the peaks themselves
            Am = A_peak[tone_i, :, :]

            # Noise Estimate - TODO: What is the num_sensors - num_sources normalization from?
            σ_c[tone_i] = real(tr( (I - Am * pinv(Am)) * Σ_sensors[tone_i, :, :] / (num_sensors - num_sources) ))
        end

        ## Store norms
        l1_errors[iter_i] = norm(γ - γ_prev, 1) / norm(γ, 1)

        if iter_i > all_cfgs.Convergence.min_iter
            if l1_errors[iter_i] < Δγ_max
                Δγ_max = l1_errors[iter_i]
                γmin = γ
                iteration_L1 = iter_i
            end

            if (l1_errors[iter_i] < all_cfgs.Convergence.error || iteration_L1 + all_cfgs.Convergence.delay <= iter_i)
                if all_cfgs.SBL.iter_flag 
                    println("Converged in iteration: ", iter_i, " with L1 error: ", l1_errors[iter_i])
                end

                # Exit the loop and set up returns
                num_iters = iter_i
                break   
            elseif iter_i == all_cfgs.Convergence.max_iter
                # Passed all iterations without converging better than our error threshold
                @warn "SBL did not converge within the maximum number of iterations specified." iteration=iter_i l1_error=l1_errors[iter_i]
                num_iters = iter_i
            elseif iter_i != all_cfgs.Convergence.max_iter && all_cfgs.SBL.iter_flag && (all_cfgs.SBL.status_report_num_iters % iter_i == 0)
                # Give a status report every status_report_num_iters iterations
                println("Status Report at iteration: ", iter_i, " with L1 error: ", l1_errors[iter_i])
            end
        end  # end if iter_i > all_cfgs.Convergence.min_iter

    end ## END OF ITERATIONS ##

    ## Posterior Distribution ##
    # x_post is the posterior unbiased mean 
    for tone_i in 1:num_tones
        A_i = @view(A[tone_i, :, :])
        x_post[tone_i, :, :] = repeat(γ, 1, num_snapshots) .* (A_i' / (σ_c[tone_i] * I + A_i * (repeat(γ, 1, num_sensors) .* A_i') ) * data[tone_i, :, :])
    end

    # Function Return
    γ = γmin

    report = Dict(
        "final_iteration"           => num_iters,
        "final_iter_noise_power"    => σ_c,
        "final_iter_γ"              => γ,
        "final_iter_x_post"         => x_post,
        "iteration_min_error"       => iteration_L1,
        "error_vector"              => l1_errors[1:num_iters]
    )

    return report
end

function gamma_update!(data::Array{Complex{U}, 3}, A::Array{Complex{U}, 3}, γ_numer::Array{U, 2}, γ_denom::Array{U, 2}, σ_c::Vector{U}, γ::Vector{U}, num_tones, num_sensors, num_snapshots) where {U <: AbstractFloat}
# Function abstracts the gamma update step for each frequency and returns the 
#  numerator and denominator matrices for each frequency
    for tone_i = 1:num_tones
        # We could extract the views here for clarity but that may be less efficient.
        # Note that the @view( . ) macro creates a view into the array without copying data
        # We only use this to extract the frequency-indexed slices of the 3D arrays
        # A_i = @view(A[tone_i, :, :])
        # γ_numer_i = @view(γ_numer[tone_i, :])
        # γ_denom_i = @view(γ_denom[tone_i, :])

        # Noise Power for this tone
        σ_ci = σ_c[tone_i]

        # Identity matrix of size num_sensors is done implicitly in I variable
        # TODO: Why do we need parens around repeat * @view(A..)'?
        denom = (σ_ci * I + A[tone_i, :, :] * (repeat(γ, 1, num_sensors) .* A[tone_i, :, :]'))
        ApΣYinv = A[tone_i, :, :]' / denom 

        # Update γ for numerator and denominator in place.
        γ_numer[tone_i, :] .= sum(abs2.( ApΣYinv * data[tone_i, :, :] ), dims=2) / num_snapshots
        γ_denom[tone_i, :] .= transpose(abs.( sum( transpose(ApΣYinv) .* A[tone_i, :, :], dims = 1) ))
    end

    # Returns nothing, this is an inplace modification of our function arguments
    nothing
end

end # module SparseBayesianLearning

