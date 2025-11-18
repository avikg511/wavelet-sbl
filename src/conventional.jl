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
    xlabel!("Time [min]")
    ylabel!("DOA [deg]")

    ylims!((-40, 40))

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