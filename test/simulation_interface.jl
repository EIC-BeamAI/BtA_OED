using Test
using Random
using LinearAlgebra
using BeamlineExperiments
using BeamTracking

println("Running: Unified BtA simulation interface"); flush(stdout)
@testset "Unified BtA simulation interface" begin
    op = default_bta_operating_point()
    beam = UncoupledBeamState()
    inputs = pack_bta_inputs(op.quad_currents, op.corrector_currents, beam)
    @test length(inputs) == length(bta_input_names()) == 34
    @test uncoupled_beam_parameters(UncoupledBeamState(inputs[21:34])) == inputs[21:34]
    @test bta_input_names()[16:end] == [
        [Symbol("corrector_current_", name) for name in CORRECTOR_NAMES]...,
        :x, :px, :y, :py, :δ0, :β_x, :α_x, :β_y, :α_y, :η_x, :η_px, :ε_x, :ε_y, :σ_δ]
    @test bta_measurements(inputs) ≈ predict_profiles(default_beam_moments())

    result = simulate_bta(beam; op,
                          rng=MersenneTwister(17))
    @test result.n_particles == 1000
    @test size(result.initial_distribution) == (1000, 6)
    @test size(result.observations.noiseless) == (4, 4)
    @test result.observations.names == PROFILE_NAMES
    @test result.observations.columns == (:mean_x, :sigma_x, :mean_y, :sigma_y)
    @test result.observations.units == :mm
    @test result.inputs == inputs
    @test size(result.jacobian) == (16, 34)
    @test result.jacobian_input_names == bta_input_names()
    @test result.jacobian_reference == :lossless_moment_observations
    @test all(isfinite, result.jacobian)
    @test result.observations.noiseless == bta_measurements(inputs)
    @test result.observations.noisy != result.observations.noiseless
    @test result.n_lost == count(!, result.final_distribution.alive)
    @test result.fractional_beam_loss == result.n_lost / 1000
    @test size(result.final_distribution.coordinates) == (1000 - result.n_lost, 6)
    @test size(result.final_distribution.all_coordinates) == (1000, 6)
    @test result.final_distribution.coordinates ==
          result.final_distribution.all_coordinates[result.final_distribution.alive, :]
    @test result.final_distribution.particle_ids == findall(result.final_distribution.alive)
    @test (result.final_distribution.first_loss_element .== 0) == result.final_distribution.alive
    @test issorted(result.observations.particle_counts; rev=true)
    @test all(result.observations.particle_counts .>= size(result.final_distribution.coordinates, 1))
    @test all(isfinite, result.final_distribution.coordinates)
    @info "Unified simulation: default 1000-particle survival" result.n_lost result.fractional_beam_loss

    # The returned tangent is exactly for vec(loss-free summaries), in input units.
    # Check all 34 columns: currents, centroids and both planes' second moments.
    steps = vcat(fill(1e-2, 20), [1e-6, 1e-7, 1e-6, 1e-7,
                                1e-6, 1e-3, 1e-4, 1e-3, 1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-8])
    for column in 1:34
        plus, minus = copy(inputs), copy(inputs)
        plus[column] += steps[column]
        minus[column] -= steps[column]
        fd = vec(bta_measurements(plus) - bta_measurements(minus)) / (2steps[column])
        @test isapprox(result.jacobian[:, column], fd; rtol=3e-5, atol=2e-5)
    end
    @test all(iszero, result.jacobian[:, [14, 15, 20]]) # Downstream of all harps.
    @test norm(result.jacobian[[5, 6, 7, 8, 13, 14, 15, 16], 1]) > 1e-4

    again = simulate_bta(beam; op,
                         rng=MersenneTwister(17), compute_jacobian=false)
    @test isnothing(again.jacobian)
    @test again.observations.noisy == result.observations.noisy
    @test again.initial_distribution == result.initial_distribution
    @test isequal(again.observations.particle_noiseless, result.observations.particle_noiseless)
    @test again.final_distribution.alive == result.final_distribution.alive

    # Changing N has no effect on the ideal noisy observation for the same seed.
    smaller = simulate_bta(beam; op,
        rng=MersenneTwister(17), nparticles=10, compute_jacobian=false, noise_std_mm=0.1)
    @test smaller.observations.noisy == result.observations.noisy
    clean = simulate_bta(beam; op,
        rng=MersenneTwister(17), nparticles=10, compute_jacobian=false, noise_std_mm=0.0)
    @test clean.observations.noiseless == clean.observations.noisy
    @test clean.initial_distribution == smaller.initial_distribution
    @test isequal(clean.observations.particle_noiseless, clean.observations.particle_noisy)

    # Use paired basis samples with exactly zero mean and identity population covariance.
    normals = vcat(sqrt(5.0) * Matrix{Float64}(I, 5, 5),
                   -sqrt(5.0) * Matrix{Float64}(I, 5, 5))
    sampled = simulate_bta(beam; op,
        nparticles=10, standard_particles=normals, compute_jacobian=false,
        noise_std_mm=(0.0, 0.05, 0.0, 0.1), apertures=false,
        observations=(:mw166, :mw006, :end), rng=MersenneTwister(1))
    X = sampled.initial_distribution
    mean = vec(sum(X; dims=1)) / 10
    centred = X .- transpose(mean)
    @test mean ≈ default_beam_moments().mean
    @test transpose(centred) * centred / 10 ≈ default_beam_moments().covariance
    @test sampled.fractional_beam_loss == 0.0
    @test sampled.observations.particle_counts == [10, 10, 10]
    @test sampled.observations.names == (:mw166, :mw006, :end)
    @test sampled.observations.noisy[:, [1, 3]] == sampled.observations.noiseless[:, [1, 3]]
    @test size(sampled.observations.noise_std_mm) == (3, 4)
    direct = track_particles(X; apertures=false)
    @test sampled.final_distribution.coordinates ≈ direct.coords.v
    mx = sum(direct.coords.v[:, 1]) / 10
    my = sum(direct.coords.v[:, 3]) / 10
    end_summary = 1000 .* [mx, norm(direct.coords.v[:, 1] .- mx) / sqrt(10),
                          my, norm(direct.coords.v[:, 3] .- my) / sqrt(10)]
    @test vec(sampled.observations.particle_noiseless[3, :]) ≈ end_summary

    # A point ensemble follows the accepted reference orbit, with zero population width.
    box = control_box(:quadrupole)
    point_beam = UncoupledBeamState(; centroid=box.centroid)
    point = simulate_bta(point_beam;
        op=box.op, nparticles=1, standard_particles=zeros(1, 5),
        compute_jacobian=false, noise_std_mm=0.0)
    @test point.fractional_beam_loss == 0.0
    @test point.observations.particle_counts == ones(Int, 4)
    @test all(iszero, point.observations.particle_noiseless[:, [2, 4]])

    named_state = (; centroid=point_beam.centroid, epsilon_x=point_beam.epsilon_x,
        beta_x=point_beam.beta_x, alpha_x=point_beam.alpha_x,
        epsilon_y=point_beam.epsilon_y, beta_y=point_beam.beta_y, alpha_y=point_beam.alpha_y)
    named_result = simulate_bta(named_state;
        op=box.op, nparticles=1, standard_particles=zeros(1, 5),
        compute_jacobian=false, noise_std_mm=zeros(4, 4))
    @test pack_bta_inputs(box.op.quad_currents, box.op.corrector_currents, named_state) == point.inputs
    @test named_result.observations.noisy == point.observations.noiseless
    @test named_result.final_distribution.coordinates == point.final_distribution.coordinates

    # Full loss must be explicit, with no frozen lost states presented as endpoint particles.
    lost = simulate_bta(beam; op,
        nparticles=10, standard_particles=fill(1e6, 10, 5),
        compute_jacobian=false, noise_std_mm=0.0)
    @test lost.fractional_beam_loss == 1.0
    @test size(lost.final_distribution.coordinates) == (0, 6)
    @test lost.observations.particle_counts == zeros(Int, 4)
    @test all(isnan, lost.observations.particle_noiseless)
    @test all(isnan, lost.observations.particle_noisy)
    @test all(lost.final_distribution.first_loss_element .> 0)
    @test lost.observations.noiseless == result.observations.noiseless

    @test_throws DimensionMismatch UncoupledBeamState(ones(9))
    @test_throws ArgumentError UncoupledBeamState(; beta_x=0.0)
    @test_throws DimensionMismatch pack_bta_inputs(ones(14), zeros(5), beam)
    @test_throws DimensionMismatch bta_measurements(ones(33))
    @test_throws ArgumentError simulate_bta(beam; op,
                                          nparticles=0)
    @test_throws ArgumentError simulate_bta(beam; op,
                                          noise_std_mm=-0.1)
    @test_throws DimensionMismatch simulate_bta(beam; op,
        nparticles=10, standard_particles=zeros(10, 6))
end
