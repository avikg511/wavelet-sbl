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
using Wavelets
using WaveletsExt

function sbl(all_cfgs::PhysicalConfigs{T, U, V}, plotting::Bool = true) where {T <: Integer, U <: AbstractFloat, V <: Complex}
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

    if !plotting
        return report
    end

    # Sample Plot of final γ values (x-axis = angles in degrees)
    angles = all_cfgs.Geometry.angles
    plt = plot(angles, report["final_iter_γ"], xlabel="DOA [deg]", ylabel="Final γ Values", title="SBL Final γ Values per Angle")
    display(plt)
    savefig(pwd() * "/presentation/angle_vs_gamma.png")
    return report
end

function wavelet_denoise_data(all_cfgs::PhysicalConfigs{T, U, V}) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    new_data = all_cfgs.Experiment.data

    # The goal is to denoise the wavelet data to simplify the SBL process.
    # I'm also curious to see how this changes the number of iterations to converge.
    wt = wavelet(WT.sym10)
    levels = 3

    # Thresholds
    thr = 0.34

    # We also need to process all the data. The data is currently (63, 878, 9)
    #   or [N sensors x T Snapshots x M tones]. We should denoise each sensor individually
    #   , then threshold, and then reconstruct our signals before passing it into SBL
    num_tones = all_cfgs.Experiment.num_tones
    num_sensors = all_cfgs.Geometry.num_sensors - length(all_cfgs.Geometry.broken_indices)

    # Create dummy signal to set things up with. We also need to process the real
    #   and imaginary data separately. Wavelets are linear though so this is fine
    xᵣ = similar(new_data[1, :, 1], U)          # ᵣ for real
    xᵢ = similar(new_data[1, :, 1], U)          # ᵢ for imaginary
    for i in 1:num_sensors
        for j in 1:num_tones
            # Wavelet denoising will be done on a per-sensor-tone basis. Let's extract
            # each sensor's data stream.
            xᵣ = real(new_data[i, :, j])
            xᵢ = imag(new_data[i, :, j])

            # Wavelet denoise for each of the streams
            xr_dwt = dwt!(xᵣ, real(new_data[i, :, j]), wt, levels)
            map!(x -> abs(x) < thr ? 0.0 : x, xr_dwt, xr_dwt)#  xᵣ, xᵣ)
            xᵣ = idwt(xr_dwt, wt, levels)

            xi_dwt = dwt!(xᵢ, imag(new_data[i, :, j]), wt, levels)
            map!(x -> abs(x) < thr ? 0.0 : x, xi_dwt, xi_dwt)
            xᵢ = idwt(xi_dwt, wt, levels)   # Temporarily store in xₜ

            # Replace original signal
            new_data[i, :, j] = xᵣ + 1im .* xᵢ
        end
    end

    # Data is now fine, let's replace it. ExperimentConfig is immutable so reconstruct
    #   a new one, hopefully we just pass in a pointer to data rather than copying
    exp_cfg = ExperimentConfig{T, U, V}(
        num_tones = all_cfgs.Experiment.num_tones,
        num_sensors = all_cfgs.Experiment.num_sensors,
        num_snapshots = all_cfgs.Experiment.num_snapshots,

        start_time = all_cfgs.Experiment.start_time,
        stop_time = all_cfgs.Experiment.stop_time,
        data = new_data,
        all_tones = all_cfgs.Experiment.all_tones
    )

    all_cfgs.Experiment = exp_cfg

    return all_cfgs
end

function sbl_denoising(all_cfgs::PhysicalConfigs{T, U, V}) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    # Denoise the data with wavelets (threshold and remove really small coefficients)
    #   and then sbl
    all_cfgs = wavelet_denoise_data(all_cfgs)

    # Now apply sbl and see results
    report = sbl(all_cfgs, false)

    # Sample Plot of final γ values
    angles = all_cfgs.Geometry.angles
    plt = plot(angles, report["final_iter_γ"], xlabel="Angle (DoA)", ylabel="Final γ Values", title="Wavelet Denoised Data, SBL Final γ Values per Angle")
    display(plt)
    savefig(pwd() * "/presentation/denoised_angle_vs_gamma.png")

    # Now wavelet processing on the final gamma with symlets after the processing is complete. The goal here is to see if any peaks were removed entirely.

    # Extract original γ values
    γ_orig = normalize(report["final_iter_γ"])
    γ_wavelet = similar(γ_orig)

    # Apply wavelet denoising to γ values
    wt = wavelet(WT.sym4)
    thr = 0.03

    # Treat γ as a signal across angles and apply DWT with thresholding
    γ_coeffs = dwt(vec(γ_orig), wt)
    γ_thresh = map(x -> abs2(x) < thr ? 0.0 : x, γ_coeffs)
    γ_wavelet .= idwt(γ_thresh, wt)

    # Handle near-zero values for plotting (avoid log(0))
    γ_wavelet = map(x -> isapprox(x, 0.0) ? 1e-9 : x, γ_wavelet)

    # Combined plot - Original and Wavelet-denoised overlaid
    angles = all_cfgs.Geometry.angles
    plt = plot(angles, γ_orig, label="Wavelet Denoised", xlabel="Angle (DoA)", ylabel="Final γ Values", title="SBL Final γ Values: Wavelet Denoised vs Wavelet Gamma Denoised", linewidth=2)
    plot!(plt, angles, γ_wavelet, label="Wavelet Denoised Wavelet Gamma", linewidth=2)
    display(plt)
    savefig(pwd() * "/presentation/denoised_gamma_with_waveletgamma.png")

    return report
end

function sbl_with_symlets(all_cfgs::PhysicalConfigs{T, U, V}) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    # Call regular sbl first, and then do wavelet processing in post
    report = sbl(all_cfgs)

    # Extract original γ values
    γ_orig = normalize(report["final_iter_γ"])
    γ_wavelet = similar(γ_orig)

    # Apply wavelet denoising to γ values
    wt = wavelet(WT.sym4)
    thr = 0.03

    # Treat γ as a signal across angles and apply DWT with thresholding
    γ_coeffs = dwt(vec(γ_orig), wt)
    γ_thresh = map(x -> abs2(x) < thr ? 0.0 : x, γ_coeffs)
    γ_wavelet .= idwt(γ_thresh, wt)

    # Handle near-zero values for plotting (avoid log(0))
    γ_wavelet = map(x -> isapprox(x, 0.0) ? 1e-9 : x, γ_wavelet)

    # Combined plot - Original and Wavelet-denoised overlaid
    angles = all_cfgs.Geometry.angles
    plt = plot(angles, γ_orig, label="Original", xlabel="Angle (DoA)", ylabel="Final γ Values", title="SBL Final γ Values: Original vs Wavelet Denoised", linewidth=2)
    plot!(plt, angles, γ_wavelet, label="Wavelet Denoised", linewidth=2)
    display(plt)
    savefig(pwd() * "/presentation/sbl_gamma_comparison.png")

    # Update report with wavelet-denoised γ
    report["final_iter_γ_wavelet"] = γ_wavelet

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
