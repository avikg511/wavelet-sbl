#=  @file conventional.jl
    @author Avik Ghosh
    @date November 18th, 2025

    Purpose: Implementation of Conventional Beamforming algorithms for Conventional Beamforming
                using simple steering vectors

    More detail about each implementation can be found in the respective functions
=#

module ConventionalBeamforming

using ..SwellConfigs
using Plots
using Wavelets

function conventional_bf(all_cfgs::PhysicalConfigs{T, U, V}) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    println("In the Conventional Beamforming Module!")

    # Conventional Beamforming Processing done for each tone
    beamforming_output = zeros(V, length(all_cfgs.Geometry.angles), all_cfgs.Experiment.num_snapshots, all_cfgs.Experiment.num_tones)
    for tone_index in 1:all_cfgs.Experiment.num_tones
        beamforming_output[:, :, tone_index] = conv_bf_per_tone(all_cfgs, tone_index)
    end

    # Sample plot
    tone_ind::T = 9
    plt = heatmap(1:all_cfgs.Experiment.num_snapshots, all_cfgs.Geometry.angles, 10*log10.(abs.(beamforming_output[:, :, tone_ind]) .^ 2))

    title!("Beamforming for Tone $tone_ind")
    xlabel!("Time [min]")
    ylabel!("DOA [deg]")

    ylims!((-40, 40))

    savefig(plt, pwd() * "/presentation/conventional_bf_op_$tone_ind.png")
    display(plt)

    #=
        Wavelet Domain Processing!
        =============================================================
        Now, what we have is beamforming output, so for M=9 tones, we have M matrices that are angles x snapshots for the expected DoA per tone.

        We won't use the above system though. We now want to work with, each snapshot, using a matrix where each row is assigned to an angle, and then the row itself contains the intensities per tone. Now, you can treat this like a signal in and of itself, where the leakage across frequency bins makes it a poor visualization.

        Applying wavelets, splitting it into high/low frequency portions and then smoothing the components with thresholding should make the overall plots cleaner when displayed.

        Though note that the last axis will be our snapshots. This can be applied 2D but first, let's just compare the effect of wavelets on smoothing/isolating the DoA plots on a per snapshot-basis.
    =#
    wt = wavelet(WT.sym4)
    thr = 0.3

    # Currently we are at [angles x snapshots x tones]
    # We want [snapshots x angles x tones]
    reorg = permutedims(beamforming_output, [2, 1, 3])
    M = 100         # Random Snapshot of interest
    old_snap = reorg[M, :, :]
    wavelet_snap = similar(old_snap)

    for i in 1:length(all_cfgs.Geometry.angles)
        # Extract the one dimensional signal across each angle, contains all tones
        sig = old_snap[i, :]

        # Take the Discrete Wavelet Transform and threshold, remove the noise-like coefficients. Can't use 0 because that shows really poor results for the dB power heatmap
        c = dwt(sig, wt)
        @show maximum(abs2, c)

        c_thresh = map(x -> abs2(x) < thr ? 0.0im : x, c)

        try
            wavelet_snap[i, :] = idwt(c_thresh, wt)
        catch
            v = all_cfgs.Geometry.angles[i]
            @show c_thresh, "Angle was $v",
            return
        end

        # Now note that we can't use any 0 coefficients because the heatmap is just pure white, rather than dark. let's fix that using isapprox
        wavelet_snap[i, :] = map(x -> isapprox(x, 0.0) ? 1e-9 : x, wavelet_snap[i, :])
    end

    # Now plotting properly for each scenario
    # Plot 1 - No Wavelet
    plt = heatmap(1:all_cfgs.Experiment.num_tones, all_cfgs.Geometry.angles, 10*log10.(abs2.(old_snap)))

    title!("Unsmoothed Tone vs Power at DoA (Snapshot $M)")
    xlabel!("ith Tone")
    ylabel!("Degree of Arrival (DoA) [deg]")

    ylims!((-40, 40))

    savefig(plt, pwd() * "/presentation/conventional_bf_unsmoothed_tonevsangle.png")
    display(plt)

    # Plot 2 - Wavelets!
    plt = heatmap(1:all_cfgs.Experiment.num_tones, all_cfgs.Geometry.angles, 10*log10.(abs2.(wavelet_snap)))

    title!("Smoothed Tone vs Power at DoA (Snapshot $M)")
    xlabel!("ith Tone")
    ylabel!("Degree of Arrival (DoA) [deg]")

    ylims!((-40, 40))

    savefig(plt, pwd() * "/presentation/conventional_bf_smoothed_tonevsangle.png")
    display(plt)
end

function conv_bf_per_tone(all_cfgs::PhysicalConfigs{T, U, V}, tone_index) where {T <: Integer, U <: AbstractFloat, V <: Complex}
    # Extract wavelength, which is just speed divided by the tone transmited, λ = c / f
    # Does sound speed change with depth? Can vertical sensor spacing be used to improve results?
    λ_transmit = all_cfgs.Sensor.sound_speed / all_cfgs.Experiment.all_tones[tone_index];
    ϕ_shift = 2 * pi * all_cfgs.Geometry.sensor_spacing / λ_transmit;
    N = all_cfgs.Geometry.num_sensors - length(all_cfgs.Geometry.broken_indices);
    sin_θs = sind.(all_cfgs.Geometry.angles);

    # Make a vector from 1 : N but remove elements in all_cfgs.Geometry.broken_indices
    # Then subtract 1 to make it zero indexed. TODO: This is potentially buggy in case
    #  broken indices is zero indexed already. Just following MATLAB implementation here.
    sensor_inds = setdiff(1:all_cfgs.Geometry.num_sensors, all_cfgs.Geometry.broken_indices) .- 1

    # Sensing Matrix
    sensing_mat = 1 / √(N) * exp.( 1im * ϕ_shift * (reshape(sensor_inds, :, 1) * reshape(sin_θs', 1, :)));

    # There is an estimated_srcs variable that is unused in the original MATLAB
    beamformer_op = zeros(V, length(sin_θs), all_cfgs.Experiment.num_snapshots)

    # Find the contribution in the sensing_mat directions per snapshot and normalize
    for snapshot_ind in 1:all_cfgs.Experiment.num_snapshots
        beamformer_op[:, snapshot_ind] = sensing_mat' * all_cfgs.Experiment.data[:, snapshot_ind, tone_index];
        beamformer_op[:, snapshot_ind] ./= maximum(abs.(beamformer_op[:, snapshot_ind]));
    end

    # Return the beamformer output for this tone
    return beamformer_op
end

end # module ConventionalBeamforming
